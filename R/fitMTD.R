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
#'   The weights are constrained to be non-negative, as in most applications
#'   of the model; Raftery's original formulation also allows negative weights
#'   provided all the transition probabilities stay in \eqn{[0, 1]}, which is
#'   not supported here. Under this constraint the likelihood is maximized by
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
#'   With \code{se = TRUE} the standard errors are those of the observed
#'   information at the maximum, by the delta method, obtained in closed form
#'   with the entries of each probability vector (the weights and each row of
#'   \eqn{Q}) summing to one. Entries that are on the boundary of the parameter
#'   space (below \eqn{10^{-4}}) are held fixed, as in the usual treatment of
#'   constrained maxima, and have \code{NA}; a warning is given if a weight is
#'   among them, since then the standard errors of the others are conditional on
#'   that weight being zero. They are asymptotic, and approximate in small
#'   samples or with sparse \eqn{Q}.
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
#' @param se Logical: compute standard errors (only for \code{model = "MTD"})?
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
#'   \item{iterations, converged}{EM iterations and convergence flag of the
#'     returned fit}
#'   \item{model}{\code{"MTD"} or \code{"MTDg"}}
#'   \item{se}{only if \code{se = TRUE}: a list with the standard errors
#'     \code{lambda} of the weights and \code{Q} of the matrix (same layout as
#'     \code{estimate}); \code{NA} for the entries that cannot be estimated, see
#'     Details}
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
fitMTD <- function(sequence, order = 2, start = NULL, nstart = 1,
                   tol = 1e-10, maxit = 10000L, model = c("MTD", "MTDg"),
                   accelerate = FALSE, se = FALSE) {
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
  if (se && model != "MTD") {
    stop("standard errors are only available for model = \"MTD\"")
  }
  maxit <- as.integer(maxit)

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
  best <- NULL
  fitsMTD <- vector("list", order)
  for (k in seq_len(order)) {
    pc <- .mtdPatterns(lagged, k)
    patterns <- pc$patterns
    counts <- pc$counts
    cols <- lagged[, seq_len(k + 1L), drop = FALSE]

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
      fit <- bestOf(lapply(starts, function(s0)
        tryCatch(runEM(FALSE, patterns, counts, s0$lambda, s0$Q),
                 error = function(e) NULL)))
    }
  }
  if (!fit$converged) {
    warning("the EM algorithm did not converge in ", maxit, " iterations")
  }

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
              model = model)
  if (se) {
    pc <- .mtdPatterns(lagged, order)
    seen <- seq_len(r) %in% (unique(as.vector(lagged[, -1L])) + 1L)
    out$se <- .mtdStdErrors(pc$patterns, pc$counts, as.numeric(lambda), Q, seen,
                            states, names(lambda))
  }
  out
}

# distinct (x_t, x_{t-1}, ..., x_{t-k}) patterns of the first k + 1 columns of
# `lagged` (0-based states) and their frequencies
.mtdPatterns <- function(lagged, k) {
  cols <- lagged[, seq_len(k + 1L), drop = FALSE]
  key <- do.call(paste, c(as.data.frame(cols), sep = "\r"))
  first <- !duplicated(key)
  list(patterns = cols[first, , drop = FALSE],
       counts = as.numeric(table(factor(key, levels = key[first]))))
}

# r x r matrix [from, to] of the weighted transitions x_{t-g} -> x_t
.lagCounts <- function(patterns, counts, g, r) {
  idx <- patterns[, g + 1L] * r + patterns[, 1L] + 1L
  v <- numeric(r * r)
  tot <- rowsum(counts, idx)
  v[as.integer(rownames(tot))] <- tot
  matrix(v, r, r, byrow = TRUE)
}

# EM for the weights only, the matrices Qs (rows are departure states) being fixed
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

# Standard errors of the MTD estimates from the observed information, by the
# delta method. In every block that must sum to one (the weights and each row
# of Q) the entries that are not on the boundary (above eps) are free, except
# the largest one, which is one minus the others; the entries at the boundary
# are held fixed and get NA. The log-likelihood is a sum of log f, with f
# bilinear in (lambda, Q), so its Hessian is available in closed form.
.mtdStdErrors <- function(patterns, counts, lambda, Q, seen, states, lagNames,
                          eps = 1e-4) {
  k <- length(lambda)
  r <- nrow(Q)
  m <- nrow(patterns)
  P <- k + r * r
  phi <- c(lambda, as.vector(Q))
  blocks <- c(list(seq_len(k)),
              lapply(which(seen), function(i) k + i + r * (seq_len(r) - 1L)))
  cols <- list()
  for (b in blocks) {
    act <- b[phi[b] > eps]
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
                Q = matrix(NA_real_, r, r, dimnames = list(states, states)))
  if (length(cols) == 0L) return(empty)
  A <- do.call(cbind, cols)
  # gradient of f_p and the cross terms of its Hessian
  G <- matrix(0, m, P)
  f <- numeric(m)
  S <- matrix(0, P, P)
  w <- counts
  for (g in seq_len(k)) {
    cell <- patterns[, g + 1L] + 1L + r * patterns[, 1L]
    q <- Q[cbind(patterns[, g + 1L] + 1L, patterns[, 1L] + 1L)]
    f <- f + lambda[g] * q
    G[, g] <- q
    G[cbind(seq_len(m), k + cell)] <- G[cbind(seq_len(m), k + cell)] + lambda[g]
  }
  if (any(f <= 0)) return(empty)
  for (g in seq_len(k)) {
    cell <- patterns[, g + 1L] + 1L + r * patterns[, 1L]
    tot <- rowsum(counts / f, cell)
    idx <- k + as.integer(rownames(tot))
    S[g, idx] <- S[g, idx] + tot
    S[idx, g] <- S[idx, g] + tot
  }
  info <- crossprod(G * (sqrt(counts) / f)) - S      # minus the Hessian
  J <- crossprod(A, info %*% A)
  covPsi <- tryCatch(solve(J), error = function(e) NULL)
  if (is.null(covPsi) || any(diag(covPsi) < 0)) {
    warning("the information matrix is singular: standard errors are not available")
    return(empty)
  }
  variance <- rowSums((A %*% covPsi) * A)
  se <- sqrt(pmax(variance, 0))
  se[rowSums(abs(A)) == 0] <- NA_real_
  if (any(lambda < eps)) {
    warning("some lag weights are on the boundary (below ", eps,
            "): their standard errors are not available (NA)")
  }
  list(lambda = stats::setNames(se[seq_len(k)], lagNames),
       Q = matrix(se[-seq_len(k)], r, r, dimnames = list(states, states)))
}
