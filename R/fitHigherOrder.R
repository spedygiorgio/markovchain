#' @title Higher order Markov Chains class
#' @name HigherOrderMarkovChain-class
#' @description The S4 class that describes \code{HigherOrderMarkovChain} objects.
#' 
#' @export
# S4 class describing a higher-order Markov chain: its states, its order and a list with
# one transition matrix per lag.
setClass("HigherOrderMarkovChain", #class name
         representation(
           states = "character", 
           order = "numeric",
           transitions = "list", 
           name = "character"
         )
#          , prototype(states = c("a","b"), byrow = TRUE, # prototypizing
#                    transitionMatrix=matrix(data = c(0,1,1,0),
#                                            nrow=2, byrow=TRUE, dimnames=list(c("a","b"), c("a","b"))),
#                    name="Unnamed Markov chain")
)

# objective function to pass to solnp: the squared distance between the
# stationary distribution X and its image sum_i lambda_i * Q_i X under the
# lambda-weighted mixture of the lag matrices (Ching et al., 2008). The value
# is divided by `scale` (its value at the starting point) because the distance
# is typically of the order of 1e-7 on real sequences, below solnp's default
# tolerance, which would otherwise stop the optimizer at the starting point.
# Objective of the least-squares ('lsq') method. For a candidate vector of weights `params`,
# the fitted distribution is sum_i params[i] * Q_i X, where Q_i X is the distribution X pushed through
# the lag-i matrix. The function returns the squared distance to the observed X.
# QX, X and scale are not arguments: they are found in the environment set in fitHigherOrder().
.fn1=function(params)
{
  QX <- get("QX")
  X <- get("X")
  scale <- get("scale")
  fitted <- 0
  for (i in 1:length(QX)) {
    fitted <- fitted + params[i] * QX[[i]]
  }
  return(sum((fitted - X)^2) / scale)
}

# The constraint 'weights sum to 1' (solnp compares this sum with eqB = 1).
# equality constraint function to pass to solnp
.eqn1=function(params){
  return(sum(params))
}

#' @name fitHigherOrder
#' @aliases seq2freqProb seq2matHigh 
#' @title Functions to fit a higher order Markov chain
#'
#' @description Given a sequence of states arising from a stationary state, it
#'   fits the underlying Markov chain distribution with higher order.
#' @usage  
#' fitHigherOrder(sequence, order = 2, method = c("mle", "lsq"))
#' seq2freqProb(sequence)
#' seq2matHigh(sequence, order)
#'
#' @param sequence The sequence of states: a character vector, or a vector
#'   coercible to character (numbers and factors are matched to the states as
#'   character strings), without missing values. A list of such vectors is
#'   taken as independent sequences, see Details.
#' @param order Markov chain order
#' @param method How the weights \eqn{\lambda} are estimated: \code{"mle"}
#'   (default) or \code{"lsq"}, see Details.
#' @return A list containing lambda, Q, and X.
#'
#' @details The fitted model expresses the distribution of the next state as
#'   the mixture \eqn{\sum_{i=1}^{k} \lambda_i Q_i x_{t-i}} of the empirical
#'   lag-\eqn{i} transition matrices \eqn{Q_i} (see \code{seq2matHigh}), with
#'   weights \eqn{\lambda_i \ge 0} summing to one. The matrices \eqn{Q_i}
#'   are the same for both methods; only the weights differ.
#'
#'   \code{method = "lsq"} (the default before version 1.2, and the only method
#'   before the argument existed) chooses \eqn{\lambda} to minimize the squared distance
#'   between the stationary distribution and its image under the mixture, as in
#'   Ching et al.; it needs the \pkg{Rsolnp} package and returns \code{NULL}
#'   with a warning if it is unavailable. This criterion is weak: each lag
#'   matrix maps the empirical distribution onto itself up to end effects,
#'   \eqn{Q_i X \approx X} with an error of order \eqn{i/n}, so the
#'   objective is nearly flat in \eqn{\lambda} and the weights it returns can
#'   be unstable from one sample to another. It remains available, to
#'   reproduce earlier results or the method of Ching et al., but it is no
#'   longer the default.
#'
#'   \code{method = "mle"} (the default) chooses \eqn{\lambda} to maximize the
#'   log-likelihood \eqn{\sum_{t=k+1}^{n} \log \sum_i \lambda_i Q_i[x_t,
#'   x_{t-i}]} of the observations a model of order \eqn{k} can predict. For
#'   fixed \eqn{Q_i} the problem is concave, so the maximum is global, and it
#'   is solved by the EM algorithm for mixture weights, without
#'   \pkg{Rsolnp}. The weights are therefore those that give the highest
#'   value of \code{\link{higherOrderLogLik}} for the same observations.
#'   The iteration stops when the log-likelihood increases by less than
#'   \eqn{10^{-10}}, or after 10000 iterations; the latter happens when a
#'   weight tends to 0 or the likelihood is nearly flat in some direction
#'   (lags carrying the same information), where EM converges slowly, and a
#'   warning then reports that the weights are approximate.
#'   They are maximum likelihood estimates \emph{conditional on} the empirical
#'   matrices \eqn{Q_i}, which are not re-estimated: this is not the maximum
#'   likelihood estimator of the mixture with free matrices, in which the
#'   weights are in general not identifiable (a common distribution can be
#'   moved from \eqn{\lambda_j Q_j} to \eqn{\lambda_i Q_i} without changing
#'   any transition probability), so the weights should not be read as the
#'   relative importance of the lags beyond this conditional sense. Note
#'   that this is not the mixture transition distribution model of Raftery
#'   (1985), in which a single matrix is shared by all lags and is estimated
#'   together with the weights; that model is fitted by \code{\link{fitMTD}}.
#'
#'   With a list of sequences the lag matrices \eqn{Q_i} and the state
#'   frequencies are computed from the counts pooled over the sequences, which
#'   are treated as independent: no pair of observations spans two sequences,
#'   and a sequence with no more than \code{order} states contributes to the
#'   frequencies only. The weights are then estimated as for a single sequence,
#'   from the observations that a model of that order can predict in every
#'   sequence.
#'
#'   A warning is given when the sequence supplies no more observations than the
#'   model has free parameters (the number used by \code{\link{higherOrderLogLik}}),
#'   since the weights are then essentially undetermined.
#'
#' @references 
#' Ching, W. K., Huang, X., Ng, M. K., & Siu, T. K. (2013). Higher-order markov 
#' chains. In Markov Chains (pp. 141-176). Springer US.
#' 
#' Ching, W. K., Ng, M. K., & Fung, E. S. (2008). Higher-order multivariate
#' Markov chains and their applications. Linear Algebra and its Applications,
#' 428(2), 492-507.
#'
#' Raftery, A. E. (1985). A model for high-order Markov chains. Journal of the
#' Royal Statistical Society, Series B, 47(3), 528-539.
#'
#' @author Giorgio Spedicato, Tae Seung Kang

#'
#' @examples
#' sequence<-c("a", "a", "b", "b", "a", "c", "b", "a", "b", "c", "a", "b",
#'             "c", "a", "b", "c", "a", "b", "a", "b")
#' # weights by maximum likelihood (the default, no Rsolnp needed)
#' fit <- fitHigherOrder(sequence, order = 2)
#' fit$lambda
#' # weights by least squares on the stationary distribution, as in Ching et al.
#' fitHigherOrder(sequence, order = 2, method = "lsq")
#' higherOrderLogLik(sequence, fit)$logLik
#'
#' @export
# Main function. The model is  P(x_t = j | past) = sum_i lambda_i * Q_i[j, x_{t-i}] : a weighted
# average (mixture) of one ordinary transition matrix per lag i = 1..order. Fitting has two parts:
# (1) the matrices Q_i are the EMPIRICAL lag-i transition matrices (counted from the data);
# (2) the weights lambda are estimated by maximum likelihood ('mle') or least squares ('lsq').
fitHigherOrder<-function(sequence, order = 2, method = c("mle", "lsq")) {
  # Step 0: choose the method (the first value, 'mle', is the default).
  method <- match.arg(method)
  # one sequence, or a list of independent sequences; as in fitMTD(), numbers
  # and factors are matched to the states as character
  # Step 1: validation. A single sequence is turned into a list of one sequence, so the code below
  # always works with lists.
  seqs <- .asSequenceList(
    sequence, 2L, "sequence must be a vector of at least two states without missing values")
  lens <- lengths(seqs)
  if (sum(lens) < 2L) {
    stop("sequence must be a vector of at least two states without missing values")
  }
  if (length(order) != 1L || is.na(order) || !is.finite(order) ||
      order < 1 || order != floor(order) || order >= max(lens)) {
    stop("order must be a positive integer smaller than the sequence length")
  }
  order <- as.integer(order)
  # Step 2: warn if there are too few observations for the number of parameters.
  nStates <- length(unique(unlist(seqs)))
  .warnFewObservations(
    nobs = sum(pmax(lens - order, 0L)),
    npar = (nStates - 1L) * (1L + order * (nStates - 1L)),
    what = paste0("a higher order Markov chain of order ", order))
  # The 'lsq' method needs the optional package Rsolnp: if it is missing we warn and return NULL.
  if (method == "lsq" && !requireNamespace("Rsolnp", quietly = TRUE)) {
    warning("package Rsolnp is required for method = \"lsq\" and is not available: ",
            "NULL is returned. Use method = \"mle\", which does not need it.",
            call. = FALSE)
    return(NULL)
  }
  # Step 3: X = frequency of each state, Q = list of the lag-1, lag-2, ... transition matrices.
  # probability of each state of the sequence(s) and lag-o transition matrices
  X <- seq2freqProb(unlist(seqs))
  Q <- lapply(seq_len(order), function(o) .seq2matHighList(seqs, o))
  # Step 4a: 'mle' -> EM algorithm for the weights (function below).
  if (method == "mle") return(.fitHigherOrderMle(seqs, Q, X, order))
  # Step 4b: 'lsq' -> constrained least squares with solnp: weights >= 0 (LB) and summing to 1.
  QX <- lapply(Q, function(Qo) Qo %*% X)
  environment(.fn1) <- environment()
  # Start from equal weights; `scale` normalises the objective so solnp does not stop immediately
  # (the distance is tiny on real data).
  params <- rep(1/order, order)
  scale <- 1
  scale0 <- .fn1(params)
  if (is.finite(scale0) && scale0 > 0) scale <- scale0
  model <- Rsolnp::solnp(params, fun=.fn1, eqfun=.eqn1, eqB=1, 
                         LB=rep(0, order), control=list(trace=0))
  list(lambda = model$pars, Q = Q, X = X)
}

# Helper: always return the data as a list of character vectors (each one an independent sequence).
# Returns `sequence` as a list of character vectors. `sequence` is either an
# atomic vector (a single sequence, which must have at least `minLength`
# states, otherwise `msgSingle` is the error) or a list of atomic vectors
# (independent sequences). No missing values are allowed.
.asSequenceList <- function(sequence, minLength, msgSingle) {
  if (is.list(sequence)) {
    ok <- length(sequence) > 0L &&
      all(vapply(sequence, function(s) is.atomic(s) && !is.null(s), logical(1)))
    if (!ok || anyNA(unlist(sequence, use.names = FALSE))) {
      stop("sequence must be a vector of states or a list of such vectors, without missing values")
    }
    return(lapply(sequence, as.character))
  }
  if (!is.atomic(sequence) || length(sequence) < minLength || anyNA(sequence)) stop(msgSingle)
  list(as.character(sequence))
}

# Pooled lag-`order` transition matrix over several sequences. Pairs (x[t-order], x[t]) are counted
# INSIDE each sequence only, then the counts of all sequences are added up and each column is
# normalised (Q[to, from] stores probabilities by column).
# Lag-`order` transition matrix Q[to, from] pooled over a list of independent
# sequences: the pairs (x_{t-order}, x_t) are counted within each sequence, so
# no pair crosses from one sequence to the next. For a single sequence it is
# seq2matHigh(); a state never observed as `from` gets a uniform column, as
# there.
.seq2matHighList <- function(seqs, order) {
  if (length(seqs) == 1L) return(seq2matHigh(seqs[[1L]], order))
  states <- names(seq2freqProb(unlist(seqs, use.names = FALSE)))
  r <- length(states)
  counts <- matrix(0, r, r)
  for (s in seqs) {
    n <- length(s)
    if (n <= order) next
    # idx = state numbers; tabulate() counts how many times each (from, to) pair occurs, in one pass.
    idx <- match(s, states)
    counts <- counts + matrix(tabulate(idx[(order + 1L):n] + r * (idx[seq_len(n - order)] - 1L),
                                       nbins = r * r), r, r)
  }
  # Columns of states never observed as 'from' have no information: give them a uniform distribution.
  colTotals <- colSums(counts)
  out <- sweep(counts, 2L, ifelse(colTotals == 0, 1, colTotals), "/")
  out[, colTotals == 0] <- 1 / r
  dimnames(out) <- list(states, states)
  out
}

# Warning helper used by the fitting functions.
# Warns when a model is fitted on no more observations than it has free
# parameters: the estimates (the weights, in particular, which then stay at
# their starting values) cannot be trusted. Nothing else changes.
.warnFewObservations <- function(nobs, npar, what) {
  if (nobs <= npar) {
    warning("only ", nobs, " observation", if (nobs != 1) "s" else "",
            " can be used to estimate the ", npar, " free parameters of ", what,
            ": the estimates are not reliable", call. = FALSE)
  }
  invisible(NULL)
}

# EM algorithm (Expectation-Maximisation) for the weights. Intuition: for each observation, how
# much is each lag responsible for the probability of what happened? 'responsibility of lag o' =
# lambda_o * q_o / (sum of all). The new weight of lag o is the average responsibility. Repeat until
# the log-likelihood stops improving; it never decreases.
# Maximum likelihood weights of the higher order model. The lag-o transition
# matrices Q_o are the empirical ones (seq2matHigh), exactly as for the least
# squares fit, and the weights maximize
#   sum_{t = order + 1}^{n} log( sum_o lambda_o * Q_o[x_t, x_{t-o}] )
# over the simplex (the sum running over the observations of every sequence a
# model of this order can predict). For fixed Q_o the objective is concave in
# lambda, so there are no local maxima, and the EM iteration for mixture
# weights increases it monotonically. The probability of every observed pair
# is positive (it is a count divided by a column total), hence the mixture
# probability of every observation is positive and the iteration is well
# defined.
#
# The likelihood depends on an observation only through its history
# (x_t, x_{t-1}, ..., x_{t-order}), and a chain with r states has at most
# r^(order + 1) of them, usually far fewer than the observations. The
# observations with the same history are therefore grouped and each distinct
# history enters the sums with its number of occurrences as weight: the EM
# iterates are the same as with one term per observation (up to rounding in
# the last digits), at a fraction of the cost on long sequences.
#
# Near a boundary of the simplex (a weight tending to 0, as when a lag adds
# nothing to the others) EM converges slowly, and maxit iterations may not
# reach tol: a warning says so, and the weights returned are those of the
# last iteration, whose log-likelihood is at most marginally below the
# maximum (EM never decreases it).
.fitHigherOrderMle <- function(seqs, Q, X, order, tol = 1e-10, maxit = 10000L) {
  states <- rownames(Q[[1L]])
  r <- length(states)
  # hist[t, ] = (x_t, x_{t-1}, ..., x_{t-order}) as state numbers, one row per
  # observation a model of this order can predict.
  hist <- do.call(rbind, lapply(seqs[lengths(seqs) > order], function(s) {
    idx <- match(s, states)
    times <- (order + 1L):length(idx)
    vapply(0:order, function(o) idx[times - o], integer(length(times)))
  }))
  if (is.null(hist)) hist <- matrix(integer(0), 0L, order + 1L)
  nobs <- nrow(hist)
  # A key identifying the history: the number sum_o (x_{t-o} - 1) r^o, exact in
  # double precision while r^(order + 1) <= 2^53, otherwise a string.
  key <- if (r^(order + 1) <= 2^53) {
    drop((hist - 1L) %*% r^(0:order))
  } else {
    do.call(paste, c(as.data.frame(hist), sep = " "))
  }
  # Distinct histories, in order of first appearance, and their counts.
  first <- which(!duplicated(key))
  counts <- tabulate(match(key, key[first]), nbins = length(first))
  hist <- hist[first, , drop = FALSE]
  # q[h, o] = probability that lag o gives to the arrival state of history h,
  # Q_o[x_t, x_{t-o}].
  q <- matrix(0, nrow = length(first), ncol = order)
  for (o in seq_len(order)) q[, o] <- Q[[o]][cbind(hist[, 1L], hist[, o + 1L])]
  # Start from equal weights.
  lambda <- rep(1 / order, order)
  logLikOld <- -Inf
  converged <- FALSE
  for (iter in seq_len(maxit)) {
    # E-step: weighted[h, o] = lambda_o * q[h, o];  mixture[h] = total probability of history h.
    weighted <- sweep(q, 2L, lambda, "*")
    mixture <- rowSums(weighted)
    logLik <- sum(counts * log(mixture))
    gain <- logLik - logLikOld
    if (gain < tol) {
      converged <- TRUE
      break
    }
    logLikOld <- logLik
    # M-step: new lambda_o = average over the observations of the responsibilities
    # weighted[h, o] / mixture[h], each history counted as many times as it occurs.
    lambda <- colSums(counts * weighted / mixture) / nobs
  }
  if (!converged) {
    warning("the EM iteration for the weights of the order ", order,
            " model stopped after ", maxit, " iterations before converging",
            " (log-likelihood still increasing by ", signif(gain, 3),
            " per iteration): a weight is probably close to 0, and the weights",
            " returned are approximate", call. = FALSE)
  }
  list(lambda = lambda, Q = Q, X = X)
}


#' Log-likelihood, deviance and information criteria of a higher order Markov chain
#'
#' @description Evaluates the log-likelihood of an empirical sequence under the
#'   higher order Markov chain returned by \code{\link{fitHigherOrder}}, and
#'   derives the deviance, AIC and BIC, so that models of different orders can
#'   be compared.
#'
#' @details The fitted model is the mixture-transition-distribution model of
#'   Raftery (1985) in the form used by Ching et al.: the probability of moving
#'   to state \eqn{x_t} given the past is
#'   \deqn{P(x_t \mid x_{t-1}, \dots, x_{t-k}) = \sum_{i=1}^{k} \lambda_i\, Q_i[x_t, x_{t-i}],}
#'   where \eqn{Q_i} is the lag-\eqn{i} transition matrix (see
#'   \code{\link{seq2matHigh}}) and \eqn{k} is the order. The log-likelihood is
#'   the sum of the logarithms of these probabilities over the observations
#'   \eqn{t = } \code{start}, \eqn{\dots}, \eqn{T}, and the deviance is
#'   \eqn{-2} times the log-likelihood.
#'
#'   Two points matter when interpreting the output. First,
#'   \code{fitHigherOrder} chooses \eqn{\lambda} by default
#'   (\code{method = "mle"}) to maximize this log-likelihood, conditionally on
#'   the empirical lag matrices, for the observations a model of that order can
#'   predict; with \code{method = "lsq"} the weights come from a least squares
#'   fit on the stationary distribution, not from maximum likelihood, so the
#'   value returned is the log-likelihood \emph{of the fitted model}, not the
#'   maximum attainable one. Second, a model of order \eqn{k} can only be evaluated from
#'   observation \eqn{k + 1} onwards; to compare orders on exactly the same data
#'   set \code{start} to \eqn{1 +} the largest order compared, otherwise the
#'   models are evaluated on different numbers of observations and neither the
#'   log-likelihood nor the information criteria are comparable.
#'
#'   The number of parameters used for AIC and BIC is the dimension of the
#'   set of transition laws the model can represent,
#'   \eqn{(r - 1)(1 + k (r - 1))}, with \eqn{r} the number of states: each of
#'   the \eqn{k} lag matrices has \eqn{r (r - 1)} free probabilities and there
#'   are \eqn{k - 1} free weights, but the weights of a mixture of lag
#'   matrices are not identifiable (a distribution common to all departure
#'   states can be moved from \eqn{\lambda_j Q_j} to \eqn{\lambda_i Q_i}
#'   without changing any transition probability), which removes
#'   \eqn{r (k - 1)} parameters from the naive count
#'   \eqn{k\, r (r - 1) + (k - 1)}: for each next state the transition
#'   probability is a sum of one term for each lag, a main-effects function of
#'   the \eqn{k} past states, with \eqn{1 + k (r - 1)} free coefficients.
#'   For \eqn{k = 1} both counts are \eqn{r (r - 1)}. For a fit returned by
#'   \code{\link{fitMTD}}, whose lags share a single matrix, it is
#'   \eqn{r (r - 1) + (k - 1)}.
#'
#' @param sequence The empirical sequence of states, a character vector or a
#'   vector coercible to character (numbers and factors are matched to the
#'   states of the fit as character strings), or a list of such vectors, taken
#'   as independent sequences: the log-likelihood is then the sum over the
#'   sequences.
#' @param fit The list returned by \code{\link{fitHigherOrder}} (or by
#'   \code{\link{fitMTD}}) for \code{sequence}. If \code{NULL}, \code{fitHigherOrder(sequence, order)} is
#'   computed.
#' @param order Order of the model to fit when \code{fit} is \code{NULL}
#'   (ignored otherwise; the order is then \code{length(fit$lambda)}).
#' @param start Index of the first observation included in the likelihood.
#'   Defaults to \code{order + 1}, the earliest observation a model of that
#'   order can predict. For a list of sequences it applies to each sequence,
#'   and the sequences with fewer than \code{start} observations are skipped.
#'
#' @return A list with components \code{logLik}, \code{deviance}, \code{AIC},
#'   \code{BIC}, \code{nobs} (number of observations entering the likelihood),
#'   \code{npar}, \code{order} and \code{start}. The log-likelihood is
#'   \code{-Inf} if an observed transition has probability zero under the model
#'   (possible only for a sequence other than the one the model was fitted on).
#'
#' @references
#' Raftery, A. E. (1985). A model for high-order Markov chains. Journal of the
#' Royal Statistical Society, Series B, 47(3), 528-539.
#'
#' Ching, W. K., Huang, X., Ng, M. K., & Siu, T. K. (2013). Higher-order markov
#' chains. In Markov Chains (pp. 141-176). Springer US.
#'
#' @seealso \code{\link{fitHigherOrder}}, \code{\link{fitMTD}},
#'   \code{\link{higherOrderPredict}}
#'
#' @examples
#' sequence <- c("a", "a", "b", "b", "a", "c", "b", "a", "b", "c", "a", "b",
#'               "c", "a", "b", "c", "a", "b", "a", "b")
#' # compare orders 1 and 2 on the same observations (start = 3)
#' fit1 <- fitHigherOrder(sequence, order = 1)
#' fit2 <- fitHigherOrder(sequence, order = 2)
#' sapply(list(order1 = fit1, order2 = fit2), function(f)
#'   unlist(higherOrderLogLik(sequence, f, start = 3)[c("logLik", "deviance", "AIC", "BIC")]))
#'
#' @export
# Log-likelihood and information criteria of a fitted higher-order model. The first `start - 1`
# observations are skipped because they have no full past; use the same `start` to compare models of
# different orders (they must be evaluated on the same observations).
higherOrderLogLik <- function(sequence, fit = NULL, order = 2, start = NULL) {
  isList <- is.list(sequence)
  # numeric or factor sequences are matched to the states as character, as in fitMTD()
  seqs <- .asSequenceList(sequence, 2L, "sequence must be a vector of states without missing values")
  # If no fitted model is given, fit one now.
  if (is.null(fit)) {
    fit <- fitHigherOrder(sequence, order)
    if (is.null(fit)) stop("package Rsolnp is required to fit the model")
  }
  if (!is.list(fit) || is.null(fit$lambda) || !is.list(fit$Q) ||
      length(fit$Q) != length(fit$lambda)) {
    stop("fit must be the list returned by fitHigherOrder()")
  }
  lambda <- as.numeric(fit$lambda)
  k <- length(lambda)
  states <- rownames(fit$Q[[1L]])
  if (is.null(states)) stop("fit$Q must have state names")
  if (anyNA(match(unlist(seqs, use.names = FALSE), states))) {
    stop("sequence contains states unknown to the fitted model")
  }
  if (is.null(start)) start <- k + 1L
  if (length(start) != 1L || !is.finite(start) || start != floor(start) ||
      start < k + 1L || (!isList && start > length(seqs[[1L]]))) {
    stop("start must be an integer between order + 1 and the sequence length")
  }
  # for a list, the sequences with no observation from `start` onwards are skipped
  seqs <- seqs[lengths(seqs) >= start]
  if (length(seqs) == 0L) {
    stop("no sequence has observations from position start onwards")
  }
  # For every observation from `start` onwards: probability of what was observed under the model,
  # the mixture sum_i lambda_i * Q_i[x_t, x_{t-i}] ; then add the logs.
  p <- unlist(lapply(seqs, function(s) {
    idx <- match(s, states)
    times <- start:length(s)
    pr <- numeric(length(times))
    for (i in seq_len(k)) {
      pr <- pr + lambda[i] * fit$Q[[i]][cbind(idx[times], idx[times - i])]
    }
    pr
  }), use.names = FALSE)
  logLik <- if (any(p <= 0)) -Inf else sum(log(p))
  # Number of parameters (see the comment below), then AIC = -2 logLik + 2 npar,
  # BIC = -2 logLik + log(n) npar.
  r <- length(states)
  # fitMTD(): one matrix shared by all lags. fitHigherOrder(): one matrix per
  # lag; the weights of such a mixture are not identifiable (a distribution
  # common to all "from" states can be moved between lambda_i Q_i and
  # lambda_j Q_j), and the set of transition laws it represents has dimension
  # (r - 1)(1 + k (r - 1)) = k r (r - 1) + (k - 1) - r (k - 1)
  npar <- if (identical(fit$model, "MTD")) r * (r - 1) + (k - 1) else (r - 1) * (1 + k * (r - 1))
  list(logLik = logLik, deviance = -2 * logLik,
       AIC = -2 * logLik + 2 * npar, BIC = -2 * logLik + log(length(p)) * npar,
       nobs = length(p), npar = npar, order = k, start = as.integer(start))
}
