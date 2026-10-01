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
#'   MTD model can have several local maxima (Berchtold, 2001), so it is
#'   advisable to use more than one starting point (\code{nstart}). The first
#'   starting point is deterministic (equal weights and the transition matrix
#'   estimated from all lags pooled); the others are drawn at random, so call
#'   \code{\link{set.seed}} beforehand for reproducible results. The fit with
#'   the highest likelihood is returned.
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
#'   A state that never occurs in a conditioning position (for instance one
#'   observed only at the end of the sequence) does not enter the likelihood;
#'   its row of \eqn{Q} is not identified and is returned as uniform.
#'
#'   AIC and BIC use the number of free parameters \eqn{r(r-1) + k - 1}, and
#'   BIC the number of observations entering the likelihood. Berchtold and
#'   Raftery (2002) do not count the elements of \eqn{Q} estimated as exactly
#'   zero; their BIC can be obtained by subtracting the number of those
#'   elements from \code{npar}.
#'
#' @param sequence An empirical sequence of states (a character vector, or a
#'   vector coercible to character), without missing values.
#' @param order Order of the model, a positive integer smaller than the length
#'   of the sequence.
#' @param start Index of the first observation entering the likelihood, at
#'   least \code{order + 1} (the default).
#' @param nstart Number of starting points of the EM algorithm (the first is
#'   deterministic, the others random).
#' @param tol Convergence tolerance on the relative change of the
#'   log-likelihood between two iterations.
#' @param maxit Maximum number of EM iterations for each starting point.
#'
#' @return A list with components
#'   \item{lambda}{the estimated lag weights, named \code{lag1}, \code{lag2}, \dots}
#'   \item{estimate}{a \code{markovchain} object (by row) with the estimated
#'     transition matrix \eqn{Q}}
#'   \item{Q}{a list of \code{order} copies of \eqn{Q} stored by column
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
#'   \item{model}{the string \code{"MTD"}}
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
#' @seealso \code{\link{fitHigherOrder}}, \code{\link{higherOrderLogLik}},
#'   \code{\link{markovchainFit}}
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
#' @export
fitMTD <- function(sequence, order = 2, start = NULL, nstart = 1,
                   tol = 1e-10, maxit = 10000L) {
  if (!is.atomic(sequence) || length(sequence) < 3L || anyNA(sequence)) {
    stop("sequence must be a vector of at least three states without missing values")
  }
  sequence <- as.character(sequence)
  n <- length(sequence)
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

  states <- sort(unique(sequence))
  r <- length(states)
  if (r < 2L) stop("sequence must contain at least two different states")
  idx <- match(sequence, states) - 1L
  times <- start:n
  # distinct patterns (x_t, x_{t-1}, ..., x_{t-order}) and their counts
  lagged <- vapply(0:order, function(g) idx[times - g], integer(length(times)))
  if (!is.matrix(lagged)) lagged <- matrix(lagged, nrow = 1L)
  key <- do.call(paste, c(as.data.frame(lagged), sep = "\r"))
  first <- !duplicated(key)
  patterns <- lagged[first, , drop = FALSE]
  counts <- as.numeric(table(factor(key, levels = key[first])))

  # deterministic start: equal weights and the transitions of all lags pooled
  pooled <- matrix(0, r, r)
  for (g in seq_len(order)) {
    pairs <- tabulate(lagged[, g + 1L] * r + lagged[, 1L] + 1L, nbins = r * r)
    pooled <- pooled + matrix(pairs, r, r, byrow = TRUE)
  }
  rowTotals <- rowSums(pooled)
  pooled[rowTotals == 0, ] <- 1
  Q0 <- pooled / rowSums(pooled)

  best <- NULL
  for (s in seq_len(nstart)) {
    if (s == 1L) {
      lambda0 <- rep(1 / order, order)
      Qs <- Q0
    } else {
      lambda0 <- stats::rexp(order)
      lambda0 <- lambda0 / sum(lambda0)
      Qs <- matrix(stats::rexp(r * r), r, r)
      Qs <- Qs / rowSums(Qs)
    }
    fit <- .mtdEM(patterns, counts, lambda0, Qs, tol, as.integer(maxit))
    if (is.null(best) || fit$logLik > best$logLik) best <- fit
  }
  if (!best$converged) {
    warning("the EM algorithm did not converge in ", maxit, " iterations")
  }

  Q <- best$Q
  # rows of states never observed in a conditioning position do not enter the
  # likelihood and are not identified: report them as uniform, whatever the start
  conditioning <- unique(as.vector(lagged[, -1L]))
  Q[setdiff(seq_len(r) - 1L, conditioning) + 1L, ] <- 1 / r
  dimnames(Q) <- list(states, states)
  lambda <- setNames(as.numeric(best$lambda), paste0("lag", seq_len(order)))
  npar <- r * (r - 1L) + order - 1L
  nobs <- length(times)
  logLik <- best$logLik
  list(lambda = lambda,
       estimate = new("markovchain", states = states, byrow = TRUE,
                      transitionMatrix = Q, name = paste0("MTD(", order, ")")),
       Q = rep(list(t(Q)), order),
       X = seq2freqProb(sequence),
       logLikelihood = logLik, npar = npar,
       AIC = -2 * logLik + 2 * npar, BIC = -2 * logLik + log(nobs) * npar,
       nobs = nobs, order = order, start = start,
       iterations = best$iterations, converged = best$converged,
       model = "MTD")
}
