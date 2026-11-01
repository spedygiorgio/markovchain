context("closestReversible(norm = \"frobenius\") and a free stationary distribution")

s4 <- letters[1:4]
P4 <- matrix(c(0.1, 0.6, 0.2, 0.1, 0.05, 0.15, 0.7, 0.1,
               0.3, 0.05, 0.25, 0.4, 0.5, 0.2, 0.1, 0.2), 4, byrow = TRUE,
             dimnames = list(s4, s4))
mc4 <- new("markovchain", states = s4, transitionMatrix = P4)
frob2 <- function(R) sum((R - P4)^2)
detailedBalance <- function(R, pi) max(abs(pi * R - t(pi * R)))

test_that("the default norm = \"pi\" is unchanged", {
  pi <- as.numeric(steadyStates(mc4))
  additive <- (P4 + sweep(t(P4), 2L, pi, "*") / pi) / 2
  res <- closestReversible(mc4)
  expect_equal(res$chain@transitionMatrix, additive, check.attributes = FALSE)
  expect_identical(res$norm, "pi")
  expect_null(res$starts)
  expect_named(res, c("chain", "stationaryDistribution", "distance", "frobeniusDistance", "norm"))
})

test_that("the fixed-pi Frobenius problem agrees with PyDTMC and is optimal", {
  skip_if_not_installed("quadprog")
  # PyDTMC 10.0.0, MarkovChain.closest_reversible(initial_distribution):
  # uniform distribution (its default) and stationary distribution of P4
  pydtmc <- list(
    uniform = c(0.11249091413878613, 0.33122651439468814, 0.22501874157716564, 0.33126382988936004,
                0.3312265143946882, 0.15001426494581777, 0.34376316083232877, 0.17499605982716523,
                0.22501874157716584, 0.343763160832329, 0.18748689979660152, 0.24373119779390365,
                0.3312638298893601, 0.17499605982716523, 0.24373119779390348, 0.25000891248957113),
    stationary = c(0.10516093350320765, 0.30548972713111627, 0.2790594253237671, 0.31028991404190903,
                   0.30630328166408904, 0.10747082417800873, 0.44670694240660752, 0.13951895175129472,
                   0.21161303853856661, 0.33784180639210659, 0.25906161356586543, 0.19148354150346136,
                   0.34134157088905437, 0.15307338607044893, 0.27778401272890763, 0.22780103031158916))
  pis <- list(uniform = rep(0.25, 4), stationary = as.numeric(steadyStates(mc4)))
  set.seed(7)
  for (k in names(pis)) {
    res <- closestReversible(mc4, pis[[k]], norm = "frobenius")
    R <- res$chain@transitionMatrix
    Rpy <- matrix(pydtmc[[k]], 4, byrow = TRUE)
    expect_equal(as.numeric(t(R)), pydtmc[[k]], tolerance = 1e-4, info = k)
    # the exact solution is at least as good as PyDTMC's iterative one
    expect_lte(frob2(R), frob2(Rpy) + 1e-12)
    expect_lt(detailedBalance(R, pis[[k]]), 1e-15)
    expect_equal(unname(rowSums(R)), rep(1, 4))
    expect_equal(res$distance, sqrt(frob2(R)))
    expect_identical(res$norm, "frobenius")
    # no random reversible chain for the same pi is closer, nor any point on
    # the segment between it and the solution (convexity)
    gaps <- vapply(1:1000, function(b) {
      Fl <- matrix(rexp(16), 4, 4)
      Fl <- (Fl + t(Fl)) / 2               # symmetric flows
      diag(Fl) <- 0
      Z <- Fl / pis[[k]]                    # Z_ij = F_ij / pi_i
      Z <- Z / max(rowSums(Z))              # keep the off-diagonal row sums <= 1
      diag(Z) <- 1 - rowSums(Z)
      c(far = frob2(Z) - frob2(R), near = frob2(R + 1e-3 * (Z - R)) - frob2(R))
    }, numeric(2))
    expect_true(all(gaps["far", ] >= 0), info = k)
    expect_true(all(gaps["near", ] >= -1e-14), info = k)
  }
})

test_that("an already reversible chain is returned unchanged", {
  skip_if_not_installed("quadprog")
  pi <- c(0.1, 0.2, 0.3, 0.4)
  Fl <- matrix(c(0, 0.02, 0.03, 0.04, 0.02, 0, 0.05, 0.06,
                 0.03, 0.05, 0, 0.1, 0.04, 0.06, 0.1, 0), 4, byrow = TRUE)
  R <- Fl / pi
  diag(R) <- 1 - rowSums(R)
  dimnames(R) <- list(s4, s4)
  rev <- new("markovchain", states = s4, transitionMatrix = R)
  expect_equal(closestReversible(rev, norm = "frobenius")$distance, 0, tolerance = 1e-12)
  set.seed(1)
  free <- closestReversible(rev, "free", norm = "frobenius", nstart = 2)
  expect_lt(free$distance, 1e-6)
})

test_that("pi need not be stationary for the Frobenius problem", {
  skip_if_not_installed("quadprog")
  pi <- c(0.4, 0.3, 0.2, 0.1)
  expect_error(closestReversible(mc4, pi), "not stationary")
  res <- closestReversible(mc4, pi, norm = "frobenius")
  expect_lt(detailedBalance(res$chain@transitionMatrix, pi), 1e-15)
  expect_equal(as.numeric(steadyStates(res$chain)), pi, tolerance = 1e-10)
  # a reducible chain is accepted when pi is given
  red <- new("markovchain", states = s4, transitionMatrix = diag(4) * 0.5 + 0.5 *
               matrix(c(0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0), 4, byrow = TRUE,
                      dimnames = list(s4, s4)))
  expect_false(is.irreducible(red))
  expect_error(closestReversible(red, norm = "frobenius"), "irreducible")
  expect_lt(closestReversible(red, rep(0.25, 4), norm = "frobenius")$distance, 1e-12)
})

test_that("a free pi is never worse than a fixed one and beats a grid search", {
  skip_if_not_installed("quadprog")
  set.seed(42)
  free <- closestReversible(mc4, "free", norm = "frobenius")
  expect_lte(free$distance, closestReversible(mc4, norm = "frobenius")$distance)
  expect_lte(free$distance, closestReversible(mc4, rep(0.25, 4), norm = "frobenius")$distance)
  expect_true(is.reversible(free$chain))
  expect_equal(as.numeric(steadyStates(free$chain)), unname(free$stationaryDistribution),
               tolerance = 1e-10)
  expect_equal(nrow(free$starts), 5L)
  expect_equal(free$distance, min(free$starts$distance))
  # three states: the best fixed-pi solution on a grid of the simplex
  s3 <- c("x", "y", "z")
  P3 <- matrix(c(0.2, 0.7, 0.1, 0.1, 0.3, 0.6, 0.6, 0.1, 0.3), 3, byrow = TRUE,
               dimnames = list(s3, s3))
  mc3 <- new("markovchain", states = s3, transitionMatrix = P3)
  grid <- seq(0.02, 0.96, by = 0.02)
  best <- Inf
  for (a in grid) for (b in grid) {
    if (a + b >= 0.99) next
    d <- closestReversible(mc3, c(a, b, 1 - a - b), norm = "frobenius")$distance
    best <- min(best, d)
  }
  free3 <- closestReversible(mc3, "free", norm = "frobenius", nstart = 3)
  expect_lte(free3$distance, best + 1e-10)
  expect_equal(nrow(closestReversible(mc3, "free", norm = "frobenius", nstart = 1)$starts), 1L)
})

test_that("column-stochastic chains are handled", {
  skip_if_not_installed("quadprog")
  col <- new("markovchain", states = s4, byrow = FALSE, transitionMatrix = t(P4))
  expect_equal(closestReversible(col, norm = "frobenius")$chain@transitionMatrix,
               t(closestReversible(mc4, norm = "frobenius")$chain@transitionMatrix))
  expect_false(closestReversible(col, norm = "frobenius")$chain@byrow)
})

test_that("invalid combinations are rejected", {
  expect_error(closestReversible(mc4, "free"), "frobenius")
  expect_error(closestReversible(mc4, "other", norm = "frobenius"), "free")
  expect_error(closestReversible(mc4, norm = "manhattan"))
  skip_if_not_installed("quadprog")
  expect_error(closestReversible(mc4, "free", norm = "frobenius", nstart = 0), "nstart")
})
