# Reference implementation in R of the block coordinate ascent of
# fitMTD(negativeWeights = TRUE), as it was before it moved to C++
# (src/mtdFreeWeights.cpp). Used only by the tests, to check that the C++
# code follows the same steps.
.mtdFreeWeightsR <- function(patterns, counts, lambda, Q, rows, tol, maxit) {
  k <- length(lambda)
  r <- nrow(Q)
  to <- patterns[, 1L] + 1L
  from <- patterns[, -1L, drop = FALSE] + 1L
  comp <- function(Q) matrix(Q[cbind(as.vector(from), rep(to, k))], ncol = k)
  C <- comp(Q)
  f <- drop(C %*% lambda)
  logLik <- function(f) if (any(f <= 0)) -Inf else sum(counts * log(f))
  L <- logLik(f)
  if (!is.finite(L) || !markovchain:::.mtdFeasible(lambda, Q, rows)) {
    stop("the starting point of the MTD fit with negative weights is not feasible")
  }

  # Largest feasible step along a direction (by bisection: halve the interval repeatedly).
  # largest feasible step in [0, upper] along a direction, the feasible set
  # being convex
  maxStep <- function(upper, feasible) {
    if (feasible(upper)) return(upper)
    lo <- 0
    hi <- upper
    for (it in 1:60) {
      mid <- (lo + hi) / 2
      if (feasible(mid)) lo <- mid else hi <- mid
    }
    lo
  }
  # Best step length along a direction, again by bisection on the derivative of the log-likelihood
  # (concave along any line, so the derivative changes sign only once).
  # maximizer of sum c log(f + delta h) on [0, upper], concave in delta, with a
  # positive derivative at 0
  lineMax <- function(h, upper) {
    d1 <- function(delta) sum(counts * h / (f + delta * h))
    if (d1(upper) >= 0) return(upper)
    lo <- 0
    hi <- upper
    for (it in 1:60) {
      mid <- (lo + hi) / 2
      if (d1(mid) > 0) lo <- mid else hi <- mid
    }
    lo
  }
  # f stays positive along f + delta h for delta below this bound
  positiveBound <- function(h) {
    neg <- h < 0
    if (!any(neg)) Inf else min(f[neg] / -h[neg]) * (1 - 1e-12)
  }

  sweeps <- 0L
  converged <- FALSE
  while (sweeps < maxit) {
    sweeps <- sweeps + 1L
    Lstart <- L
    # weights
    w <- counts / f
    grad <- colSums(C * w)
    a <- which.max(grad)
    b <- which.min(grad)
    if (grad[a] - grad[b] > 1e-12 * abs(L)) {
      h <- C[, a] - C[, b]
      upper <- min(positiveBound(h), 10)
      upper <- maxStep(upper, function(d) {
        l <- lambda
        l[a] <- l[a] + d
        l[b] <- l[b] - d
        markovchain:::.mtdFeasible(l, Q, rows)
      })
      if (upper > 0) {
        delta <- lineMax(h, upper)
        lambda[a] <- lambda[a] + delta
        lambda[b] <- lambda[b] - delta
        f <- f + delta * h
      }
    }
    # rows of Q
    for (i in rows) {
      w <- counts / f
      # derivative of f_p with respect to Q[i, j]: sum of the weights of the
      # lags at which the departure state is i, if the arrival state is j
      lamI <- drop((from == i) %*% lambda)
      grad <- vapply(seq_len(r), function(j) sum((w * lamI)[to == j]), numeric(1))
      a <- which.max(grad)
      canDecrease <- which(Q[i, ] > 0)
      b <- canDecrease[which.min(grad[canDecrease])]
      if (length(b) == 0L || grad[a] - grad[b] <= 1e-12 * abs(L)) next
      h <- lamI * ((to == a) - (to == b))
      upper <- min(Q[i, b], positiveBound(h))
      upper <- maxStep(upper, function(d) {
        Qn <- Q
        Qn[i, a] <- Qn[i, a] + d
        Qn[i, b] <- max(Qn[i, b] - d, 0)
        markovchain:::.mtdFeasible(lambda, Qn, rows)
      })
      if (upper <= 0) next
      delta <- lineMax(h, upper)
      Q[i, a] <- Q[i, a] + delta
      Q[i, b] <- max(Q[i, b] - delta, 0)
      f <- f + delta * h
    }
    C <- comp(Q)
    f <- drop(C %*% lambda)
    L <- logLik(f)
    if (L - Lstart <= tol * (abs(L) + tol)) {
      converged <- TRUE
      break
    }
  }
  list(lambda = lambda, Q = Q, logLik = L, iterations = sweeps, converged = converged)
}
