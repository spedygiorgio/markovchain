library(testthat)
library(markovchain)

states3 <- c("a", "b", "c")

directedCycle <- new("markovchain", states = states3,
  transitionMatrix = matrix(c(0, 1, 0,
                              0, 0, 1,
                              1, 0, 0), byrow = TRUE, nrow = 3,
                            dimnames = list(states3, states3)))

triangleWalk <- new("markovchain", states = states3,
  transitionMatrix = matrix(c(0, 0.5, 0.5,
                              0.5, 0, 0.5,
                              0.5, 0.5, 0), byrow = TRUE, nrow = 3,
                            dimnames = list(states3, states3)))

# The pi-weighted Hilbert-Schmidt norm the construction minimizes.
piNorm <- function(A, pi) sqrt(sum(sweep(pi * A^2, 2L, pi, FUN = "/")))

randomChain <- function(n, seed) {
  set.seed(seed)
  M <- matrix(runif(n * n), n, n)
  M <- M / rowSums(M)
  nm <- letters[seq_len(n)]
  dimnames(M) <- list(nm, nm)
  new("markovchain", states = nm, transitionMatrix = M)
}

test_that("the approximation is a valid, reversible chain with the same stationary distribution", {
  for (n in c(2, 3, 5, 8)) {
    mc <- randomChain(n, seed = 100 + n)
    result <- closestReversible(mc)

    expect_true(validObject(result$chain), info = paste("n =", n))
    expect_true(is.reversible(result$chain), info = paste("n =", n))
    expect_equal(as.numeric(steadyStates(result$chain)),
                 as.numeric(result$stationaryDistribution),
                 tolerance = 1e-8, info = paste("n =", n))
  }
})

test_that("the closest reversible chain to a directed cycle is the undirected walk", {
  result <- closestReversible(directedCycle)
  expect_false(is.reversible(directedCycle))
  expect_equal(result$chain@transitionMatrix,
               triangleWalk@transitionMatrix, tolerance = 1e-12)
})

test_that("an already reversible chain is returned unchanged, at distance zero", {
  result <- closestReversible(triangleWalk)
  expect_equal(result$chain@transitionMatrix, triangleWalk@transitionMatrix,
               tolerance = 1e-12)
  expect_equal(result$distance, 0, tolerance = 1e-10)
  expect_equal(result$frobeniusDistance, 0, tolerance = 1e-10)
})

test_that("the result matches the additive reversibilization computed by hand", {
  mc <- randomChain(5, seed = 7)
  P <- mc@transitionMatrix
  result <- closestReversible(mc)
  pi <- as.numeric(result$stationaryDistribution)

  Pstar <- matrix(0, 5, 5)
  for (i in 1:5) for (j in 1:5) Pstar[i, j] <- pi[j] * P[j, i] / pi[i]
  expected <- (P + Pstar) / 2

  expect_equal(unname(result$chain@transitionMatrix), unname(expected),
               tolerance = 1e-12)
})

test_that("no feasible competitor is closer in the pi-weighted norm", {
  # Feasible competitors are S = R + Delta/pi with Delta symmetric and with
  # zero row sums: then pi_i*S_ij stays symmetric (reversible) and the rows
  # still sum to one.
  mc <- randomChain(5, seed = 7)
  result <- closestReversible(mc)
  P <- mc@transitionMatrix
  R <- result$chain@transitionMatrix
  pi <- as.numeric(result$stationaryDistribution)
  best <- piNorm(P - R, pi)

  set.seed(99)
  checked <- 0
  for (rep in 1:300) {
    B <- matrix(rnorm(25, sd = 0.02), 5, 5)
    B <- (B + t(B)) / 2
    diag(B) <- 0
    Delta <- B - diag(rowSums(B))
    S <- R + Delta / pi
    if (any(S < 0)) next
    checked <- checked + 1
    expect_gte(piNorm(P - S, pi), best - 1e-12)
  }
  expect_gt(checked, 20) # the check is only meaningful if competitors were tried
})

test_that("the residual is orthogonal to the reversible subspace (Pythagoras)", {
  mc <- randomChain(5, seed = 7)
  result <- closestReversible(mc)
  P <- mc@transitionMatrix
  R <- result$chain@transitionMatrix
  pi <- as.numeric(result$stationaryDistribution)

  set.seed(123)
  B <- matrix(rnorm(25, sd = 0.02), 5, 5)
  B <- (B + t(B)) / 2
  diag(B) <- 0
  S <- R + (B - diag(rowSums(B))) / pi

  expect_equal(piNorm(P - S, pi)^2,
               piNorm(P - R, pi)^2 + piNorm(R - S, pi)^2, tolerance = 1e-10)
})

test_that("the plain Frobenius distance is NOT what gets minimized", {
  # Documented caveat: moving along a feasible direction can reduce the plain
  # Frobenius distance while increasing the pi-weighted one that is minimized.
  mc <- randomChain(4, seed = 11)
  result <- closestReversible(mc)
  P <- mc@transitionMatrix
  R <- result$chain@transitionMatrix
  pi <- as.numeric(result$stationaryDistribution)
  g <- P - R

  bestDirection <- NULL
  bestInner <- 0
  for (i in 1:3) for (j in (i + 1):4) {
    E <- matrix(0, 4, 4)
    E[i, j] <- E[j, i] <- 1
    E[i, i] <- E[i, i] - 1
    E[j, j] <- E[j, j] - 1
    D <- E / pi
    inner <- sum(g * D)
    if (abs(inner) > abs(bestInner)) {
      bestInner <- inner
      bestDirection <- D
    }
  }
  # A non-zero Frobenius gradient in a feasible direction is exactly why R is
  # not the Frobenius minimizer.
  expect_gt(abs(bestInner), 1e-6)

  step <- bestInner / sum(bestDirection^2) / 4
  S <- R + step * bestDirection
  expect_true(all(S >= 0))

  flow <- pi * S
  expect_equal(flow, t(flow), tolerance = 1e-12)         # still reversible
  expect_lt(sqrt(sum((P - S)^2)), sqrt(sum(g^2)))        # closer in Frobenius
  expect_gt(piNorm(P - S, pi), piNorm(g, pi))            # but further in ||.||_pi
})

test_that("a user-supplied stationary distribution is accepted when it is stationary", {
  mc <- randomChain(4, seed = 21)
  pi <- as.numeric(steadyStates(mc))
  fromDefault <- closestReversible(mc)
  fromSupplied <- closestReversible(mc, stationaryDistribution = pi)
  expect_equal(fromSupplied$chain@transitionMatrix,
               fromDefault$chain@transitionMatrix, tolerance = 1e-10)

  # An unnormalized but proportional vector is normalized internally.
  fromScaled <- closestReversible(mc, stationaryDistribution = pi * 17)
  expect_equal(fromScaled$chain@transitionMatrix,
               fromDefault$chain@transitionMatrix, tolerance = 1e-10)
})

test_that("a non-stationary or invalid stationary distribution is rejected", {
  mc <- randomChain(4, seed = 21)
  expect_error(closestReversible(mc, stationaryDistribution = c(0.25, 0.25, 0.25, 0.25)),
               "not stationary")
  expect_error(closestReversible(mc, stationaryDistribution = c(0.5, 0.5, 0, 0)),
               "strictly positive")
  expect_error(closestReversible(mc, stationaryDistribution = c(0.5, 0.5)),
               "one entry per state")
  expect_error(closestReversible(mc, tolerance = -1), "tolerance")
})

test_that("column-stochastic storage is preserved and gives the transposed result", {
  mc <- randomChain(4, seed = 33)
  P <- mc@transitionMatrix
  mcCol <- new("markovchain", states = states(mc), byrow = FALSE,
               transitionMatrix = t(P))

  byRow <- closestReversible(mc)
  byCol <- closestReversible(mcCol)

  expect_false(byCol$chain@byrow)
  expect_true(validObject(byCol$chain))
  expect_equal(byCol$chain@transitionMatrix, t(byRow$chain@transitionMatrix),
               tolerance = 1e-12)
  expect_equal(byCol$distance, byRow$distance, tolerance = 1e-12)
})

test_that("closestReversible rejects reducible chains", {
  reducible <- new("markovchain", states = c("a", "b"),
                   transitionMatrix = matrix(c(1, 0, 0.5, 0.5),
                                              byrow = TRUE, nrow = 2,
                                              dimnames = list(c("a", "b"), c("a", "b"))))
  expect_error(closestReversible(reducible), "irreducible")
})

test_that("the one-state chain is its own closest reversible chain", {
  one_state <- new("markovchain", states = "a",
                   transitionMatrix = matrix(1, 1, 1, dimnames = list("a", "a")))
  result <- closestReversible(one_state)
  expect_equal(unname(result$chain@transitionMatrix), matrix(1, 1, 1))
  expect_equal(result$distance, 0, tolerance = 1e-12)
})
