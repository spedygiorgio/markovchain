#' Fit a mixture transition distribution (MTD) model
#'
#' @description Estimates by maximum likelihood the mixture transition
#'   distribution model of Raftery (1985) for a high-order Markov chain: the
#'   probability of the next state is a weighted mixture of the
#'   contributions of the last \code{order} states, all governed by the same
#'   transition matrix.
#'
#' @details For an order \eqn{k} the model is
#'   \deqn{P(X_t = j \mid X_{t-1} = i_1, \dots, X_{t-k} = i_k) =
#'     \sum_{g=1}^{k} \lambda_g\, q_{i_g j},}
#'   where \eqn{Q = (q_{ij})} is an \eqn{r \times r} transition matrix (rows
#'   are departure states) and the lag weights \eqn{\lambda_g} sum to one.
#'   It needs \eqn{r(r-1) + k - 1} parameters instead of the
#'   \eqn{r^k (r-1)} of a fully parameterized Markov chain of order \eqn{k},
#'   which makes high orders practicable (Raftery, 1985; Berchtold and
#'   Raftery, 2002). This is different from the model fitted by
#'   \code{\link{fitHigherOrder}}, which mixes a different empirical matrix
#'   for each lag and does not estimate them jointly with the weights.
#'
#'   By default the weights are constrained to be non-negative, as in most
#'   applications of the model (see below for negative weights). Under this
#'   constraint the likelihood is maximized by
#'   the EM algorithm of Lebre and Bourguignon (2008), in which the latent
#'   variable is the lag that generated each observation; it is implemented in
#'   C++. Every iteration increases the likelihood, but the likelihood of the
#'   MTD model can have several local maxima (Berchtold, 2001). The models of
#'   order \eqn{1, \dots,} \code{order} are therefore fitted in turn on the
#'   same observations, and order \eqn{k} is started from equal weights with
#'   the transition matrix of all lags pooled, and from the fit of order
#'   \eqn{k - 1} extended with a zero and with a small positive weight for the
#'   new lag. The first of these extensions has the likelihood of order
#'   \eqn{k - 1}, so the likelihood returned never decreases with the order,
#'   as it must for nested models. \code{nstart - 1} further random starting
#'   points can be added for the requested order; they are drawn at random, so
#'   call \code{\link{set.seed}} beforehand for reproducible results. The fit
#'   with the highest likelihood is returned.
#'
#'   The likelihood is conditional on the observations before \code{start}:
#'   it is the product of the transition probabilities of \eqn{x_t} for
#'   \eqn{t = } \code{start}, \eqn{\dots}, \eqn{n}. The default,
#'   \code{order + 1}, uses every observation an order-\code{order} model can
#'   predict. To compare models of different orders (or with the Markov chains
#'   of \code{\link{fitHigherOrder}} via \code{\link{higherOrderLogLik}}) use
#'   the same \code{start} for all, for instance \eqn{1 +} the largest order;
#'   Berchtold and Raftery (2002) condition on the first 14 observations,
#'   i.e. \code{start = 15}.
#'
#'   With a list of sequences the likelihood is the product over the
#'   sequences: each contributes the observations from \code{start} to its end,
#'   the earlier ones being only used as conditioning states, no lag spans two
#'   sequences, and a sequence with fewer than \code{start} states is not used
#'   (apart from the state frequencies \code{X}). \code{order} must be smaller
#'   than the length of the longest sequence.
#'
#'   A state that never occurs in a conditioning position (for instance one
#'   observed only at the end of the sequence) does not enter the likelihood;
#'   its row of \eqn{Q} is not identified and is returned as uniform.
#'
#'   A warning is given when the number of observations entering the
#'   likelihood does not exceed the number of free parameters.
#'
#'   With \code{model = "MTDg"} (Berchtold, 2001; Lebre and Bourguignon, 2008)
#'   the probability of the next state is \eqn{\sum_g \lambda_g\, q^{(g)}_{i_g j}}
#'   with a transition matrix \eqn{Q_g} for each lag. It has
#'   \eqn{(r-1)(1 + k(r-1))} free parameters, the number of the models of
#'   \code{\link{fitHigherOrder}}, of which it is the maximum likelihood
#'   version: there the matrices are the empirical ones and only the weights are
#'   fitted. The EM algorithm is the same; order \eqn{k} is started from the
#'   MTD fit of the same order, from the empirical matrices with the best
#'   weights for them, and from the fit of order \eqn{k - 1}, so its likelihood
#'   is never below those of \code{model = "MTD"} and of the
#'   \code{method = "mle"} fit of \code{\link{fitHigherOrder}}, nor of
#'   lower orders. The row of \eqn{Q_g} of a state never seen at lag \eqn{g} is
#'   returned as uniform.
#'
#'   The weights and matrices of the MTDg model are not identified: since the
#'   probabilities only depend on the products \eqn{A_g = \lambda_g Q_g}, adding
#'   \eqn{\epsilon} times the same probability vector to every row of \eqn{A_1}
#'   and subtracting it from every row of \eqn{A_2} changes \eqn{\lambda_1},
#'   \eqn{\lambda_2}, \eqn{Q_1} and \eqn{Q_2} but no transition probability.
#'   The maximum is therefore attained on a set of dimension up to
#'   \eqn{(k - 1) r}, of which the values returned are one point; only the
#'   transition probabilities (and the likelihood) are determined, and
#'   \code{npar} counts the identified parameters. The lag weights of an MTDg
#'   fit should not be interpreted as the importance of the lags.
#'
#'   With \code{negativeWeights = TRUE} (model \code{"MTD"} only) the weights may
#'   be negative, as in Raftery's (1985) original formulation, provided every
#'   transition probability \eqn{\sum_g \lambda_g q_{i_g j}} lies in
#'   \eqn{[0, 1]} for all the combinations of departure states \eqn{i_g}: for
#'   each arrival state \eqn{j}, \eqn{\min_i q_{ij}} times the sum of the positive
#'   weights must not be below \eqn{\max_i q_{ij}} times the sum of the absolute
#'   negative weights (Raftery, 1985). The EM algorithm does not apply, and the
#'   likelihood is maximized by block coordinate ascent in the spirit of
#'   Berchtold (2001): in the weights and in each row of \eqn{Q} in turn, mass is
#'   moved from the entry with the smallest to the one with the largest partial
#'   derivative, by the step that maximizes the likelihood along that direction
#'   within the feasible region. Order \eqn{k} starts from the fit with
#'   non-negative weights of the same order and from the fit of order
#'   \eqn{k - 1} with a zero weight for the new lag, so its likelihood is never
#'   below either of them; for the requested order it also starts with each lag
#'   in turn given a negative weight. The likelihood has several local maxima,
#'   and this is slower than the EM algorithm. \code{iterations} then counts the
#'   sweeps of the ascent at the requested order. Unidentified rows of \eqn{Q} are
#'   returned as uniform, or as the mean of the other rows if a uniform row
#'   would violate the condition.
#'
#'   With \code{se = TRUE} the standard errors are those of the observed
#'   information at the maximum, by the delta method, obtained in closed form
#'   with the entries of each probability vector (the weights and each row of
#'   \eqn{Q}) summing to one. Entries that are on the boundary of the parameter
#'   space (below \eqn{10^{-4}}) are held fixed, as in the usual treatment of
#'   constrained maxima, and have \code{NA}; a warning is given if a weight is
#'   among them, since then the standard errors of the others are conditional on
#'   that weight being zero (with negative weights allowed, if Raftery's
#'   condition is binding). Standard errors are also given for the transition
#'   probabilities of every history observed in the data. For the MTDg model
#'   only these are available, the weights and matrices not being identified;
#'   they are computed with a generalized inverse of the information, which
#'   gives the same variance for any identified quantity. All are asymptotic,
#'   and approximate in small samples or with sparse \eqn{Q}.
#'
#'   AIC and BIC use the number of free parameters \eqn{r(r-1) + k - 1}, and
#'   BIC the number of observations entering the likelihood. Berchtold and
#'   Raftery (2002) do not count the elements of \eqn{Q} estimated as exactly
#'   zero; their BIC can be obtained by subtracting the number of those
#'   elements from \code{npar}.
#'
#' @param sequence An empirical sequence of states (a character vector, or a
#'   vector coercible to character), without missing values, or a list of such
#'   vectors, taken as independent sequences (see Details).
#' @param order Order of the model, a positive integer smaller than the length
#'   of the sequence.
#' @param start Index of the first observation entering the likelihood, at
#'   least \code{order + 1} (the default).
#' @param nstart \code{nstart - 1} is the number of random starting points of
#'   the EM algorithm added to the deterministic ones (see Details).
#' @param tol Convergence tolerance on the relative change of the
#'   log-likelihood between two iterations.
#' @param maxit Maximum number of EM iterations for each starting point (EM
#'   steps, if \code{accelerate = TRUE}).
#' @param model \code{"MTD"} (default), the model with a single matrix, or
#'   \code{"MTDg"}, the general model with a different matrix for each lag.
#' @param accelerate Logical: speed up the EM algorithm with the SQUAREM
#'   scheme of Varadhan and Roland (2008)? The default \code{FALSE} keeps the
#'   plain algorithm. Convergence is much faster when the likelihood is flat,
#'   and the same maximum is reached.
#' @param se Logical: compute standard errors (see Details)?
#' @param negativeWeights Logical: allow negative lag weights, as in Raftery
#'   (1985)? Only for \code{model = "MTD"}; the default \code{FALSE} keeps them
#'   non-negative.
#'
#' @return A list with components
#'   \item{lambda}{the estimated lag weights, named \code{lag1}, \code{lag2}, \dots}
#'   \item{estimate}{a \code{markovchain} object (by row) with the estimated
#'     transition matrix \eqn{Q} (for \code{model = "MTDg"}, a list with the
#'     \code{order} objects of the matrices \eqn{Q_1, \dots, Q_k})}
#'   \item{Q}{a list of \code{order} copies of \eqn{Q} (the \eqn{Q_g} for
#'     \code{"MTDg"}) stored by column
#'     (\code{Q[[g]][to, from]}), the layout returned by
#'     \code{\link{fitHigherOrder}}, so that the fit can be passed to
#'     \code{\link{higherOrderLogLik}}}
#'   \item{X}{the relative frequencies of the states in \code{sequence}}
#'   \item{logLikelihood, npar, AIC, BIC, nobs}{maximized log-likelihood, number
#'     of free parameters, information criteria and number of observations
#'     entering the likelihood}
#'   \item{order, start}{as used in the fit}
#'   \item{iterations, converged}{EM iterations (sweeps of the coordinate
#'     ascent with \code{negativeWeights = TRUE}) and convergence flag of the
#'     returned fit}
#'   \item{model}{\code{"MTD"} or \code{"MTDg"}}
#'   \item{negativeWeights}{as used in the fit}
#'   \item{se}{only if \code{se = TRUE}: a list with the standard errors
#'     \code{lambda} of the weights, \code{Q} of the matrix (same layout as
#'     \code{estimate}; \code{NA} for the entries that cannot be estimated and,
#'     for MTDg, for all the weights and matrices) and \code{probabilities} of
#'     the transition probabilities, a matrix with one row per history observed
#'     in the data (named by its states, oldest first, separated by \code{"-"},
#'     as passed to \code{\link{higherOrderPredict}}) and one column per next
#'     state; see Details}
#'
#' @references
#' Raftery, A. E. (1985). A model for high-order Markov chains. Journal of the
#' Royal Statistical Society, Series B, 47(3), 528-539.
#'
#' Berchtold, A. (2001). Estimation in the mixture transition distribution
#' model. Journal of Time Series Analysis, 22(4), 379-397.
#'
#' Berchtold, A. and Raftery, A. E. (2002). The mixture transition distribution
#' model for high-order Markov chains and non-Gaussian time series.
#' Statistical Science, 17(3), 328-356.
#'
#' Lebre, S. and Bourguignon, P.-Y. (2008). An EM algorithm for estimation in
#' the mixture transition distribution model. Journal of Statistical
#' Computation and Simulation, 78(1), 1-15.
#'
#' Varadhan, R. and Roland, C. (2008). Simple and globally convergent methods
#' for accelerating the convergence of any EM algorithm. Scandinavian Journal of
#' Statistics, 35(2), 335-353.
#'
#' @seealso \code{\link{fitHigherOrder}}, \code{\link{higherOrderLogLik}},
#'   \code{\link{higherOrderPredict}}, \code{\link{markovchainFit}}
#'
#' @examples
#' # hourly wind directions at Koeberg (Berchtold and Raftery, 2002)
#' wind <- read.csv(system.file("extdata", "koeberg_wind.csv",
#'                              package = "markovchain"))$state
#' fit <- fitMTD(wind, order = 2, start = 15)
#' fit$lambda
#' fit$estimate
#' c(logLik = fit$logLikelihood, BIC = fit$BIC)
#'
#' # several starting points guard against local maxima
#' set.seed(1)
#' fitMTD(wind, order = 3, start = 15, nstart = 5)$logLikelihood
#'
#' # standard errors, and the model with a matrix for each lag
#' fit2 <- fitMTD(wind, order = 2, start = 15, se = TRUE)
#' fit2$se$lambda
#' fitMTD(wind, order = 2, start = 15, model = "MTDg")$logLikelihood
#' @export
# Overview of the MTD (Mixture Transition Distribution) fit. Model:
#   P(x_t = j | x_{t-1} = i1, ..., x_{t-k} = ik) = sum_g lambda_g * Q[i_g, j]
# One SINGLE matrix Q (shared by all lags) and k weights lambda that sum to 1: much fewer parameters than
# a full order-k chain (r^k (r-1)). The 'MTDg' variant uses a different matrix for each lag.
# Parameters are estimated by maximum likelihood with the EM algorithm (the heavy loops are in C++,
# see src/fitMTD.cpp). The function proceeds in the following steps:
#   1. validate all the arguments;
#   2. build the table of lagged observations;
#   3. fit orders 1, 2, ..., order in turn (each order is started from the previous fit);
#   4. optionally refit with negative weights (MTD) or with one matrix per lag (MTDg);
#   5. assemble the result: weights, matrices, logLik, AIC/BIC and, if asked, standard errors.
fitMTD <- function(sequence, order = 2, start = NULL, nstart = 1,
                   tol = 1e-10, maxit = 10000L, model = c("MTD", "MTDg"),
                   accelerate = FALSE, se = FALSE, negativeWeights = FALSE) {
  # Step 1: validation (match.arg picks one of the allowed values of `model`).
  model <- match.arg(model)
  seqs <- .asSequenceList(
    sequence, 3L, "sequence must be a vector of at least three states without missing values")
  n <- max(lengths(seqs))
  if (n < 3L) {
    stop("sequence must be a vector of at least three states without missing values")
  }
  if (length(order) != 1L || is.na(order) || !is.finite(order) ||
      order < 1 || order != floor(order) || order >= n) {
    stop("order must be a positive integer smaller than the sequence length")
  }
  order <- as.integer(order)
  if (is.null(start)) start <- order + 1L
  if (length(start) != 1L || is.na(start) || !is.finite(start) ||
      start != floor(start) || start < order + 1L || start > n) {
    stop("start must be an integer between order + 1 and the sequence length")
  }
  start <- as.integer(start)
  if (length(nstart) != 1L || is.na(nstart) || nstart < 1 || nstart != floor(nstart)) {
    stop("nstart must be a positive integer")
  }
  if (length(tol) != 1L || is.na(tol) || tol <= 0) stop("tol must be positive")
  if (length(maxit) != 1L || is.na(maxit) || maxit < 1) stop("maxit must be a positive integer")
  if (length(accelerate) != 1L || is.na(accelerate) || !is.logical(accelerate)) {
    stop("accelerate must be TRUE or FALSE")
  }
  if (length(se) != 1L || is.na(se) || !is.logical(se)) stop("se must be TRUE or FALSE")
  if (length(negativeWeights) != 1L || is.na(negativeWeights) || !is.logical(negativeWeights)) {
    stop("negativeWeights must be TRUE or FALSE")
  }
  if (negativeWeights && model != "MTD") {
    stop("negative weights are only available for model = \"MTD\"")
  }
  maxit <- as.integer(maxit)

  # Step 2: states sorted alphabetically; r = number of states.
  states <- sort(unique(unlist(seqs, use.names = FALSE)))
  r <- length(states)
  if (r < 2L) stop("sequence must contain at least two different states")
  # lagged[, g + 1] is x_{t-g} for the observations t = start..n entering the
  # likelihood; the sequences of a list are independent, so the observations
  # of every sequence with at least `start` states are stacked and no lag
  # spans two sequences
  lagged <- do.call(rbind, lapply(seqs[lengths(seqs) >= start], function(s) {
    idx <- match(s, states) - 1L
    times <- start:length(s)
    out <- vapply(0:order, function(g) idx[times - g], integer(length(times)))
    if (!is.matrix(out)) out <- matrix(out, nrow = 1L)
    out
  }))
  times <- seq_len(nrow(lagged))

  # runEM(): runs one EM from a starting point (plain C++ version, or accelerated with SQUAREM).
  # one run of the EM algorithm from the given starting point, plain or SQUAREM
  # accelerated; Q is a matrix if shared (MTD), otherwise a list of matrices
  runEM <- function(shared, patterns, counts, lambda, Q) {
    if (!accelerate) {
      if (shared) return(.mtdEM(patterns, counts, lambda, Q, tol, maxit))
      return(.mtdgEM(patterns, counts, lambda, Q, tol, maxit))
    }
    out <- .accelerateEM(patterns, counts, lambda, if (shared) list(Q) else Q,
                         shared, tol, maxit)
    if (shared) out$Q <- out$Q[[1L]]
    out
  }
  # bestOf(): of several fits started from different points, keep the one with the highest
  # log-likelihood (the likelihood can have several local maxima).
  bestOf <- function(fits) {
    fits <- Filter(Negate(is.null), fits)
    if (length(fits) == 0L) stop("the EM algorithm could not be started")
    fits[[which.max(vapply(fits, function(f) f$logLik, numeric(1)))]]
  }

  # The models of order 1, ..., order are fitted in turn on the same
  # observations. Besides the deterministic start, order k also starts from the
  # fit of order k - 1 with a zero weight for the new lag, a fixed point of the
  # EM algorithm with the same likelihood, so that the likelihood returned can
  # never decrease with the order, and with a small positive weight for it.
  # Step 3: loop over the orders 1..order.
  best <- NULL
  fitsMTD <- vector("list", order)
  for (k in seq_len(order)) {
    # Instead of the n observations we store only the DISTINCT patterns (x_t, x_{t-1}, ..., x_{t-k})
    # and how many times each occurs: the likelihood only depends on these counts, so the EM is faster.
    pc <- .mtdPatterns(lagged, k)
    patterns <- pc$patterns
    counts <- pc$counts
    cols <- lagged[, seq_len(k + 1L), drop = FALSE]

    # Starting point 1: equal weights and one matrix obtained by pooling the transitions of all lags.
    # deterministic start: equal weights and the transitions of all lags pooled
    pooled <- matrix(0, r, r)
    for (g in seq_len(k)) pooled <- pooled + .lagCounts(patterns, counts, g, r)
    pooled[rowSums(pooled) == 0, ] <- 1
    starts <- list(list(lambda = rep(1 / k, k), Q = pooled / rowSums(pooled)))
    if (k > 1L) {
      starts <- c(starts,
                  list(list(lambda = c(best$lambda, 0), Q = best$Q),
                       list(lambda = c(best$lambda * (1 - 1 / k), 1 / k), Q = best$Q)))
    }
    # Optional extra random starts (nstart > 1): random weights and matrix.
    if (k == order && nstart > 1L) {
      for (s in seq_len(nstart - 1L)) {
        lambda0 <- stats::rexp(k)
        Qs <- matrix(stats::rexp(r * r), r, r)
        starts[[length(starts) + 1L]] <- list(lambda = lambda0 / sum(lambda0),
                                              Q = Qs / rowSums(Qs))
      }
    }
    best <- bestOf(lapply(starts, function(s0)
      runEM(TRUE, patterns, counts, s0$lambda, s0$Q)))
    fitsMTD[[k]] <- best
  }
  fit <- best

  # Step 4a: negative weights allowed (Raftery 1985): see .mtdFreeWeights below.
  if (negativeWeights) {
    # Raftery's model with weights of any sign. Order k starts from the fit with
    # non-negative weights of the same order and from the fit of order k - 1
    # with a zero weight for the new lag (so that the likelihood never decreases
    # with the order, nor falls below the non-negative fit), and the likelihood
    # is maximized by block coordinate ascent under Raftery's condition.
    prev <- fitsMTD[[1L]]
    for (k in seq_len(order)[-1L]) {
      pc <- .mtdPatterns(lagged, k)
      seenRows <- sort(unique(as.vector(lagged[, 1L + seq_len(k)]))) + 1L
      starts <- list(fitsMTD[[k]], list(lambda = c(prev$lambda, 0), Q = prev$Q))
      if (k == order) {
        # the likelihood has several local maxima: for the requested order also
        # start with each lag in turn given a negative weight, the matrix being
        # shrunk towards its column means so that the negative weight is allowed
        l0 <- as.numeric(fitsMTD[[k]]$lambda)
        Q0 <- fitsMTD[[k]]$Q
        Qmean <- matrix(colMeans(Q0[seenRows, , drop = FALSE]), r, r, byrow = TRUE)
        for (g in seq_len(k)) for (v in c(-0.05, -0.15, -0.3)) for (sh in c(0, 0.3, 0.6)) {
          l <- l0
          l[g] <- v
          l[-g] <- if (sum(l[-g]) > 0) l[-g] * (1 - v) / sum(l[-g]) else (1 - v) / (k - 1)
          Qs <- (1 - sh) * Q0 + sh * Qmean
          if (.mtdFeasible(l, Qs, seenRows)) starts[[length(starts) + 1L]] <- list(lambda = l, Q = Qs)
        }
      }
      prev <- bestOf(lapply(starts, function(s0)
        tryCatch(.mtdFreeWeights(pc$patterns, pc$counts, as.numeric(s0$lambda), s0$Q,
                                 seenRows, tol, maxit),
                 error = function(e) NULL)))
    }
    if (order > 1L) fit <- prev
  }

  # Step 4b: MTDg model, one transition matrix per lag. Same scheme as above.
  if (model == "MTDg") {
    # Same scheme for the model with a matrix for each lag. Order k starts from
    # the MTD fit of order k (so that the likelihood is never below it), from the
    # empirical matrix of every lag with the best weights for them (the
    # fitHigherOrder() mle model) and with equal weights, from the fit of order
    # k - 1 with a zero and a small weight for the new lag, and from random points.
    fit <- NULL
    for (k in seq_len(order)) {
      pc <- .mtdPatterns(lagged, k)
      patterns <- pc$patterns
      counts <- pc$counts
      Qemp <- lapply(seq_len(k), function(g) {
        N <- .lagCounts(patterns, counts, g, r)
        N[rowSums(N) == 0, ] <- 1
        N / rowSums(N)
      })
      starts <- list(
        list(lambda = .weightsEM(patterns, counts, Qemp, rep(1 / k, k)), Q = Qemp),
        list(lambda = rep(1 / k, k), Q = Qemp),
        list(lambda = fitsMTD[[k]]$lambda, Q = rep(list(fitsMTD[[k]]$Q), k)))
      if (k > 1L) {
        newQ <- c(fit$Q, list(Qemp[[k]]))
        starts <- c(starts,
                    list(list(lambda = c(fit$lambda, 0), Q = newQ),
                         list(lambda = c(fit$lambda * (1 - 1 / k), 1 / k), Q = newQ)))
      }
      if (k == order && nstart > 1L) {
        for (s in seq_len(nstart - 1L)) {
          lambda0 <- stats::rexp(k)
          Qs <- lapply(seq_len(k), function(g) {
            M <- matrix(stats::rexp(r * r), r, r)
            M / rowSums(M)
          })
          starts[[length(starts) + 1L]] <- list(lambda = lambda0 / sum(lambda0), Q = Qs)
        }
      }
      # Last, the empirical matrices with the exact maximum likelihood weights
      # for them (the fitHigherOrder() mle weights, computed by Newton's
      # method), so that the likelihood is never below that model even where
      # the EM for the weights above stops short of their maximum. It comes
      # after the other starts, so that it is chosen only when it gives a
      # strictly higher likelihood.
      if (k > 1L) {
        starts[[length(starts) + 1L]] <- list(lambda = .weightsNewton(patterns, counts, Qemp),
                                              Q = Qemp)
      }
      fit <- bestOf(lapply(starts, function(s0)
        tryCatch(runEM(FALSE, patterns, counts, s0$lambda, s0$Q),
                 error = function(e) NULL)))
    }
  }
  # Report a warning if the EM did not converge within `maxit` steps.
  if (!fit$converged) {
    warning("the EM algorithm did not converge in ", maxit, " iterations")
  }

  # Step 5: assemble the output. The rows of Q for states never seen as a conditioning state do not
  # affect the likelihood, so they cannot be estimated: they are reported as uniform.
  nobs <- length(times)
  # rows of states never observed in a conditioning position (at the lag, for
  # MTDg) do not enter the likelihood and are not identified: report them as
  # uniform, whatever the start
  uniformRows <- function(Q, g) {
    seen <- unique(if (g == 0L) as.vector(lagged[, -1L]) else lagged[, g + 1L])
    Q[setdiff(seq_len(r) - 1L, seen) + 1L, ] <- 1 / r
    dimnames(Q) <- list(states, states)
    Q
  }
  lambda <- stats::setNames(as.numeric(fit$lambda), paste0("lag", seq_len(order)))
  X <- seq2freqProb(unlist(seqs, use.names = FALSE))
  logLik <- fit$logLik
  if (model == "MTD") {
    Q <- uniformRows(fit$Q, 0L)
    if (negativeWeights && !.mtdFeasible(lambda, Q, seq_len(r))) {
      # with negative weights a uniform row may give negative probabilities:
      # the unidentified rows are set to the mean of the others instead, which
      # keeps every probability within the range of the identified ones
      seenRows <- sort(unique(as.vector(lagged[, -1L]))) + 1L
      unseen <- setdiff(seq_len(r), seenRows)
      Q[unseen, ] <- matrix(colMeans(Q[seenRows, , drop = FALSE]), length(unseen), r,
                            byrow = TRUE)
    }
    # Parameters: r(r-1) for the matrix Q + (k-1) free weights.
    npar <- r * (r - 1L) + order - 1L
    estimate <- new("markovchain", states = states, byrow = TRUE,
                    transitionMatrix = Q, name = paste0("MTD(", order, ")"))
    Qcol <- rep(list(t(Q)), order)
  } else {
    Qg <- lapply(seq_len(order), function(g) uniformRows(fit$Q[[g]], g))
    npar <- (r - 1L) * (1L + order * (r - 1L))
    estimate <- lapply(seq_len(order), function(g)
      new("markovchain", states = states, byrow = TRUE, transitionMatrix = Qg[[g]],
          name = paste0("MTDg(", order, "), lag ", g)))
    Qcol <- lapply(Qg, t)
  }
  .warnFewObservations(nobs, npar, paste0("an ", model, " model of order ", order))
  out <- list(lambda = lambda, estimate = estimate, Q = Qcol, X = X,
              logLikelihood = logLik, npar = npar,
              AIC = -2 * logLik + 2 * npar, BIC = -2 * logLik + log(nobs) * npar,
              nobs = nobs, order = order, start = start,
              iterations = fit$iterations, converged = fit$converged,
              model = model, negativeWeights = negativeWeights)
  # Optional standard errors: from the observed information matrix (see .mtdStdErrors).
  if (se) {
    pc <- .mtdPatterns(lagged, order)
    if (model == "MTD") {
      seen <- list(seq_len(r) %in% (unique(as.vector(lagged[, -1L])) + 1L))
      Qs <- list(Q)
    } else {
      seen <- lapply(seq_len(order), function(g) seq_len(r) %in% (unique(lagged[, g + 1L]) + 1L))
      Qs <- Qg
    }
    out$se <- .mtdStdErrors(pc$patterns, pc$counts, as.numeric(lambda), Qs, seen,
                            states, names(lambda), freeWeights = negativeWeights)
    if (model == "MTD") out$se$Q <- out$se$Q[[1L]]
  }
  out
}

# Helper: from the table of lagged observations keeps the distinct patterns of the first k lags,
# and counts how many times each one appears.
# distinct (x_t, x_{t-1}, ..., x_{t-k}) patterns of the first k + 1 columns of
# `lagged` (0-based states) and their frequencies
.mtdPatterns <- function(lagged, k) {
  cols <- lagged[, seq_len(k + 1L), drop = FALSE]
  key <- do.call(paste, c(as.data.frame(cols), sep = "\r"))
  first <- !duplicated(key)
  list(patterns = cols[first, , drop = FALSE],
       counts = as.numeric(table(factor(key, levels = key[first]))))
}

# Helper: weighted count of the transitions x_{t-g} -> x_t as an r x r matrix [from, to].
# r x r matrix [from, to] of the weighted transitions x_{t-g} -> x_t
.lagCounts <- function(patterns, counts, g, r) {
  idx <- patterns[, g + 1L] * r + patterns[, 1L] + 1L
  v <- numeric(r * r)
  tot <- rowsum(counts, idx)
  v[as.integer(rownames(tot))] <- tot
  matrix(v, r, r, byrow = TRUE)
}

# Helper: EM that estimates only the weights, with the matrices fixed (used as starting values).
# EM for the weights only, the matrices Qs (rows are departure states) being fixed
# The same maximum likelihood weights for fixed matrices, exactly, by the
# active-set Newton method of fitHigherOrder() (R/fitHigherOrder.R).
.weightsNewton <- function(patterns, counts, Qs) {
  k <- length(Qs)
  comp <- matrix(vapply(seq_len(k), function(g)
    Qs[[g]][cbind(patterns[, g + 1L] + 1L, patterns[, 1L] + 1L)],
    numeric(nrow(patterns))), ncol = k)
  as.numeric(.mixtureWeightsNewton(comp, counts))
}

.weightsEM <- function(patterns, counts, Qs, lambda, iters = 2000L, tol = 1e-13) {
  k <- length(Qs)
  comp <- matrix(vapply(seq_len(k), function(g)
    Qs[[g]][cbind(patterns[, g + 1L] + 1L, patterns[, 1L] + 1L)],
    numeric(nrow(patterns))), ncol = k)
  total <- sum(counts)
  for (i in seq_len(iters)) {
    w <- sweep(comp, 2L, lambda, "*")
    new <- colSums(counts * w / rowSums(w)) / total
    done <- max(abs(new - lambda)) < tol
    lambda <- new
    if (done) break
  }
  lambda
}

# SQUAREM: accelerates the EM. Idea: two plain EM steps show the 'direction' in which the parameters
# are moving; instead of crawling along it the algorithm jumps ahead by an extrapolated amount
# (alpha), projects the result back to valid probabilities, and falls back (backtracking) to a
# shorter jump whenever the likelihood would decrease.
# SQUAREM acceleration (Varadhan and Roland, 2008, scheme S3) of the EM
# iteration of the MTD (shared = TRUE, Qs a list with one matrix) or MTDg model
# (Qs a list of k matrices). A cycle needs two EM steps, an extrapolation that
# is projected on the parameter space, and one more EM step to evaluate it; the
# extrapolation is backtracked towards the plain double EM step until it does
# not decrease the likelihood, so every cycle is monotone. `maxit` bounds the
# number of EM steps.
.accelerateEM <- function(patterns, counts, lambda, Qs, shared, tol, maxit) {
  k <- length(lambda)
  nQ <- length(Qs)
  r <- nrow(Qs[[1L]])
  rr <- r * r
  pack <- function(lambda, Qs) c(lambda, unlist(lapply(Qs, as.vector), use.names = FALSE))
  unpack <- function(theta) {
    list(lambda = theta[seq_len(k)],
         Qs = lapply(seq_len(nQ), function(b) matrix(theta[k + (b - 1L) * rr + seq_len(rr)], r, r)))
  }
  step <- function(theta) {
    u <- unpack(theta)
    s <- if (shared) .mtdStep(patterns, counts, u$lambda, u$Qs[[1L]])
         else .mtdgStep(patterns, counts, u$lambda, u$Qs)
    list(theta = pack(s$lambda, if (shared) list(s$Q) else s$Q), logLik = s$logLik)
  }
  project <- function(theta) {
    u <- unpack(theta)
    l <- pmax(u$lambda, 0)
    l <- if (sum(l) > 0) l / sum(l) else rep(1 / k, k)
    Qs <- lapply(u$Qs, function(M) {
      M <- pmax(M, 0)
      s <- rowSums(M)
      M[s == 0, ] <- 1
      M / rowSums(M)
    })
    pack(l, Qs)
  }
  theta <- pack(lambda, Qs)
  s0 <- step(theta)
  evals <- 1L
  if (!is.finite(s0$logLik)) stop("an observed pattern has probability zero at the starting values")
  L0 <- s0$logLik
  converged <- FALSE
  while (evals < maxit) {
    th1 <- s0$theta
    s1 <- step(th1)
    evals <- evals + 1L
    L1 <- s1$logLik
    if (L1 - L0 <= tol * (abs(L1) + tol)) {
      theta <- th1
      L0 <- L1
      converged <- TRUE
      break
    }
    th2 <- s1$theta
    res <- th1 - theta
    vv <- th2 - th1 - res
    nv <- sqrt(sum(vv^2))
    alpha <- if (nv > 0) -max(1, sqrt(sum(res^2)) / nv) else -1
    repeat {
      cand <- if (alpha == -1) th2 else project(theta - 2 * alpha * res + alpha^2 * vv)
      sc <- step(cand)
      evals <- evals + 1L
      if (alpha == -1 || (is.finite(sc$logLik) && sc$logLik >= L1)) break
      alpha <- (alpha - 1) / 2
      if (alpha > -1.05) alpha <- -1
    }
    theta <- cand
    s0 <- sc
    L0 <- sc$logLik
  }
  u <- unpack(theta)
  list(lambda = u$lambda, Q = u$Qs, logLik = L0, iterations = evals, converged = converged)
}

# Standard errors via the delta method:
#  (1) the free parameters are those that are not on the boundary (zero); in each block summing to one
#      one entry is determined by the others;
#  (2) the information matrix J is minus the Hessian of the log-likelihood, available in closed form;
#  (3) the covariance of the parameters is the (generalised) inverse of J;
#  (4) the variance of any function of the parameters (as the transition probabilities) follows from
#      its gradient D: Var = D Cov D'.
# Standard errors of the MTD (one matrix in Qs) or MTDg (one matrix per lag)
# estimates from the observed information, by the delta method. In every block
# that must sum to one (the weights and each row of each matrix) the entries
# that are not on the boundary (above eps) are free, except the largest one,
# which is one minus the others; the entries at the boundary are held fixed and
# get NA. With freeWeights (negative weights allowed) no weight is on the
# boundary. The log-likelihood is a sum of log f, with f bilinear in the weights
# and the matrices, so its Hessian is available in closed form.
.mtdStdErrors <- function(patterns, counts, lambda, Qs, seen, states, lagNames,
                          freeWeights = FALSE, eps = 1e-4) {
  k <- length(lambda)
  B <- length(Qs)
  r <- nrow(Qs[[1L]])
  rr <- r * r
  m <- nrow(patterns)
  P <- k + B * rr
  phi <- c(lambda, unlist(lapply(Qs, as.vector), use.names = FALSE))
  block <- function(g) if (B == 1L) 1L else g
  rowBlocks <- unlist(lapply(seq_len(B), function(b)
    lapply(which(seen[[b]]), function(i) k + (b - 1L) * rr + i + r * (seq_len(r) - 1L))),
    recursive = FALSE)
  cols <- list()
  for (bl in c(list(seq_len(k)), rowBlocks)) {
    isWeights <- identical(bl, seq_len(k))
    act <- if (isWeights && freeWeights) bl else bl[phi[bl] > eps]
    if (length(act) < 2L) next
    ref <- act[which.max(phi[act])]
    for (o in setdiff(act, ref)) {
      a <- numeric(P)
      a[o] <- 1
      a[ref] <- -1
      cols[[length(cols) + 1L]] <- a
    }
  }
  empty <- list(lambda = stats::setNames(rep(NA_real_, k), lagNames),
                Q = lapply(seq_len(B), function(b)
                  matrix(NA_real_, r, r, dimnames = list(states, states))))
  if (length(cols) == 0L) return(empty)
  A <- do.call(cbind, cols)
  # gradient of f_p and the cross terms of its Hessian
  G <- matrix(0, m, P)
  f <- numeric(m)
  cells <- lapply(seq_len(k), function(g)
    k + (block(g) - 1L) * rr + patterns[, g + 1L] + 1L + r * patterns[, 1L])
  for (g in seq_len(k)) {
    q <- phi[cells[[g]]]
    f <- f + lambda[g] * q
    G[, g] <- q
    G[cbind(seq_len(m), cells[[g]])] <- G[cbind(seq_len(m), cells[[g]])] + lambda[g]
  }
  if (any(f <= 0)) return(empty)
  S <- matrix(0, P, P)
  for (g in seq_len(k)) {
    tot <- rowsum(counts / f, cells[[g]])
    idx <- as.integer(rownames(tot))
    S[g, idx] <- S[g, idx] + tot
    S[idx, g] <- S[idx, g] + tot
  }
  info <- crossprod(G * (sqrt(counts) / f)) - S      # minus the Hessian
  J <- crossprod(A, info %*% A)
  J <- (J + t(J)) / 2
  # In the MTDg model the weights and the matrices are not identified (mass can
  # be moved between the lags without changing any transition probability), so
  # J is singular; the transition probabilities are identified, and their
  # variance is the same with any generalized inverse of J. For the MTD model J
  # is invertible unless the estimate is degenerate.
  eig <- eigen(J, symmetric = TRUE)
  top <- max(abs(eig$values))
  # in the MTDg model the flat directions are curved in these coordinates, so
  # away from the exact maximum their curvature is small but not zero
  flat <- if (B == 1L) 1e-6 else 1e-4
  if (any(eig$values < -flat * top)) {
    warning("the observed information is not positive semi-definite: ",
            "standard errors are not available")
    return(empty)
  }
  # directions along which the likelihood is flat (up to the convergence of the
  # fit) are dropped
  keep <- eig$values > (if (B == 1L) 1e-10 else flat) * top
  singular <- !all(keep)
  covPsi <- eig$vectors[, keep, drop = FALSE] %*%
    (t(eig$vectors[, keep, drop = FALSE]) / eig$values[keep])
  if (singular && B == 1L) {
    warning("the information matrix is singular: standard errors are not available")
    return(empty)
  }

  # transition probabilities of every context observed in the data
  contexts <- unique(patterns[, -1L, drop = FALSE])
  nc <- nrow(contexts)
  probSE <- matrix(NA_real_, nc, r)
  for (j in seq_len(r)) {
    D <- matrix(0, nc, P)
    for (g in seq_len(k)) {
      cell <- k + (block(g) - 1L) * rr + contexts[, g] + 1L + r * (j - 1L)
      D[, g] <- phi[cell]
      D[cbind(seq_len(nc), cell)] <- D[cbind(seq_len(nc), cell)] + lambda[g]
    }
    DA <- D %*% A
    probSE[, j] <- sqrt(pmax(rowSums((DA %*% covPsi) * DA), 0))
  }
  # the history of a context, oldest first, as in higherOrderPredict()
  dimnames(probSE) <- list(
    apply(contexts[, k:1, drop = FALSE] + 1L, 1L, function(h) paste(states[h], collapse = "-")),
    states)
  probSE <- probSE[order(rownames(probSE)), , drop = FALSE]

  if (B > 1L) {
    return(list(lambda = empty$lambda, Q = empty$Q, probabilities = probSE))
  }
  variance <- rowSums((A %*% covPsi) * A)
  se <- sqrt(pmax(variance, 0))
  se[rowSums(abs(A)) == 0] <- NA_real_
  if (!freeWeights && any(lambda < eps)) {
    warning("some lag weights are on the boundary (below ", eps,
            "): their standard errors are not available (NA)")
  }
  if (freeWeights && .mtdFeasibilitySlack(lambda, Qs[[1L]], which(seen[[1L]])) < 1e-6) {
    warning("the estimate lies on the boundary of Raftery's condition (a transition ",
            "probability is zero for some lags): the standard errors are unreliable")
  }
  list(lambda = stats::setNames(se[seq_len(k)], lagNames),
       Q = list(matrix(se[k + seq_len(rr)], r, r, dimnames = list(states, states))),
       probabilities = probSE)
}

# Raftery's condition: with negative weights a probability sum_g lambda_g q_g may become negative.
# This function measures the 'slack' of that constraint (>= 0 means all probabilities are valid).
# Raftery's condition for the MTD model with weights of any sign: every
# probability sum_g lambda_g q[i_g, j], over all the combinations of the
# departure states i_g in `rows`, must be non-negative (they then sum to one, so
# none exceeds one). The smallest of them for arrival state j is
# min_i q[i, j] * (sum of the positive weights) - max_i q[i, j] * (sum of the
# absolute negative weights); the slack is the smallest over j.
.mtdFeasibilitySlack <- function(lambda, Q, rows) {
  Qr <- Q[rows, , drop = FALSE]
  pos <- sum(pmax(lambda, 0))
  neg <- sum(pmax(-lambda, 0))
  min(apply(Qr, 2L, min) * pos - apply(Qr, 2L, max) * neg)
}

.mtdFeasible <- function(lambda, Q, rows, tol = 1e-12) {
  all(Q >= 0) && .mtdFeasibilitySlack(lambda, Q, rows) >= -tol
}

# Maximum likelihood with weights of any sign: block coordinate ascent. At each move a small amount
# of mass is shifted between two entries of a block that sums to one, in the direction in which the
# likelihood increases most, with the step length chosen by a one-dimensional search.
# Maximum likelihood for the MTD model with weights of any sign, by block
# coordinate ascent in the spirit of Berchtold (2001): in each block that sums
# to one (the weights, then each identified row of Q) mass moves from the entry
# with the smallest to the one with the largest partial derivative, by the
# step that maximizes the likelihood along that direction (a concave function
# of the step) within the feasible interval (entries of Q non-negative,
# Raftery's condition, positive probabilities of the observed patterns).
# Every move increases the likelihood; the sweeps stop when one increases it
# by less than tol relatively, or after maxit sweeps.
.mtdFreeWeights <- function(patterns, counts, lambda, Q, rows, tol, maxit) {
  # The ascent itself is in C++ (src/mtdFreeWeights.cpp), step by step the
  # algorithm described above; here the arguments are put in the form it
  # expects and the names are restored.
  storage.mode(patterns) <- "integer"
  storage.mode(Q) <- "double"
  out <- .mtdFreeWeightsRcpp(patterns, as.numeric(counts), as.numeric(lambda), Q,
                             as.integer(rows), as.numeric(tol), as.integer(maxit))
  dimnames(out$Q) <- dimnames(Q)
  names(out$lambda) <- names(lambda)
  out
}
