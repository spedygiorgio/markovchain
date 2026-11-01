# The top-down heuristic of aggregateStates() scores candidate merges through the flows between
# macro-states (.topDownMergeScores()) and only ranks near ties with the full computation. These
# tests check that this gives exactly the result of the original algorithm, which evaluated every
# candidate in full, and that the scores differ from the divergences by a constant.

# The original top-down loop (markovchain 1.2.1), kept here as the reference.
.referenceTopDown <- function(P, pi, k) {
  n <- nrow(P)
  Q <- P
  eta <- diag(n)
  for (i in seq_len(n - k)) {
    qPi <- if (i == 1L) pi else markovchain:::.calculateInvariant(Q)
    pairs <- markovchain:::.buildPairCandidates(Q, qPi)
    divergences <- numeric(length(pairs))
    candidates <- vector("list", length(pairs))
    for (j in seq_along(pairs)) {
      res <- markovchain:::.calculateQTopDown(P, pi, eta, pairs[[j]][1L], pairs[[j]][2L])
      candidates[[j]] <- res
      divergences[j] <- markovchain:::.klDivergenceAggregate(P, pi, res$phi, res$Q)
    }
    best <- which.min(divergences)
    Q <- candidates[[best]]$Q
    eta <- candidates[[best]]$phi
  }
  list(Q = Q, phi = eta)
}

.topdownTestChains <- function() {
  set.seed(123)
  out <- list()
  dense <- function(n) { P <- matrix(stats::rexp(n * n)^2, n); P / rowSums(P) }
  for (n in c(5, 8, 12, 16)) out[[length(out) + 1L]] <- dense(n)
  # symmetric (reversible with uniform stationary distribution): many exact ties
  A <- matrix(stats::rexp(100), 10); A <- A + t(A)
  out[[length(out) + 1L]] <- A / rowSums(A)
  # lazy random walk on a cycle: degenerate spectrum and exactly tied merges
  n <- 9; Cy <- matrix(0, n, n)
  for (i in 1:n) { Cy[i, i] <- 0.5; Cy[i, (i %% n) + 1] <- 0.25; Cy[i, ((i - 2) %% n) + 1] <- 0.25 }
  out[[length(out) + 1L]] <- Cy
  # exactly lumpable chain with uniform split inside the blocks
  g <- rep(1:3, length.out = 9)
  B <- matrix(c(.6, .3, .1, .2, .5, .3, .3, .3, .4), 3, byrow = TRUE)
  L <- matrix(0, 9, 9)
  for (i in 1:9) for (b in 1:3) { idx <- which(g == b); L[i, idx] <- B[g[i], b] / length(idx) }
  out[[length(out) + 1L]] <- L
  out
}

test_that("top-down gives exactly the result of the original algorithm", {
  for (P in .topdownTestChains()) {
    n <- nrow(P)
    mc <- new("markovchain", transitionMatrix = P, states = paste0("s", seq_len(n)))
    pi <- as.numeric(steadyStates(mc))
    pi <- pi / sum(pi)
    for (k in unique(c(2L, n %/% 2L, n - 1L))) {
      expect_identical(markovchain:::.aggregateSpectralTopDown(P, pi, k),
                       .referenceTopDown(P, pi, k), info = paste("n =", n, "k =", k))
    }
  }
})

test_that("merge scores differ from the divergences by a constant", {
  set.seed(9)
  n <- 10
  P <- matrix(stats::rexp(n * n), n); P <- P / rowSums(P)
  P[sample(n * n, 20)] <- 0; P <- P / rowSums(P)          # some zero transitions
  mc <- new("markovchain", transitionMatrix = P, states = paste0("s", seq_len(n)))
  pi <- as.numeric(steadyStates(mc)); pi <- pi / sum(pi)
  # a partition into 6 macro-states, and its flows
  grp <- c(1, 1, 2, 3, 3, 3, 4, 5, 6, 6)
  eta <- outer(grp, 1:6, `==`) * 1
  W <- t(eta) %*% (pi * P) %*% eta
  Pi <- as.numeric(t(eta) %*% pi)
  pairs <- list(c(1L, 2L), c(3L, 5L), c(4L, 6L), c(2L, 6L))
  scores <- markovchain:::.topDownMergeScores(W, Pi, pairs)
  exact <- vapply(pairs, function(pr) {
    res <- markovchain:::.calculateQTopDown(P, pi, eta, pr[1L], pr[2L])
    markovchain:::.klDivergenceAggregate(P, pi, res$phi, res$Q)
  }, numeric(1))
  expect_equal(scores - scores[1], exact - exact[1], tolerance = 1e-10)
})
