library(testthat)
library(markovchain)

states2 <- c("a", "b")
P2 <- matrix(c(0.7, 0.3, 0.1, 0.9),
             byrow = TRUE, nrow = 2,
             dimnames = list(states2, states2))
mc2 <- new("markovchain", states = states2, transitionMatrix = P2)

states3 <- c("a", "b", "c")

## ---- is.reversible -------------------------------------------------------

test_that("a 2-state irreducible chain is always reversible", {
  expect_true(is.reversible(mc2))

  # A very asymmetric 2-state chain is still reversible: with 2 states,
  # detailed balance is just the stationarity equation restated.
  asym <- new("markovchain", states = states2,
              transitionMatrix = matrix(c(0.99, 0.01, 0.02, 0.98),
                                        byrow = TRUE, nrow = 2,
                                        dimnames = list(states2, states2)))
  expect_true(is.reversible(asym))
})

test_that("an undirected random walk on a triangle is reversible", {
  triangle <- new("markovchain", states = states3,
                  transitionMatrix = matrix(c(0, 0.5, 0.5,
                                              0.5, 0, 0.5,
                                              0.5, 0.5, 0), byrow = TRUE, nrow = 3,
                                            dimnames = list(states3, states3)))
  expect_true(is.reversible(triangle))
})

test_that("a directed cycle is not reversible", {
  cycle3 <- new("markovchain", states = states3,
               transitionMatrix = matrix(c(0, 1, 0,
                                           0, 0, 1,
                                           1, 0, 0), byrow = TRUE, nrow = 3,
                                         dimnames = list(states3, states3)))
  expect_false(is.reversible(cycle3))
})

test_that("is.reversible matches a direct detailed-balance check", {
  # An asymmetric, non-trivial 3-state chain.
  P3 <- matrix(c(0.5, 0.3, 0.2,
                 0.2, 0.6, 0.2,
                 0.1, 0.1, 0.8), byrow = TRUE, nrow = 3,
               dimnames = list(states3, states3))
  mc3 <- new("markovchain", states = states3, transitionMatrix = P3)
  pi <- as.numeric(steadyStates(mc3))
  flow <- outer(pi, rep(1, 3)) * P3
  expect_equal(is.reversible(mc3), isTRUE(all.equal(flow, t(flow), tolerance = 1e-8)))
})

test_that("is.reversible supports column-stochastic storage", {
  mc2_col <- new("markovchain", states = states2, byrow = FALSE,
                 transitionMatrix = t(P2))
  expect_equal(is.reversible(mc2_col), is.reversible(mc2))
})

test_that("is.reversible is true for the one-state chain", {
  one_state <- new("markovchain", states = "a",
                   transitionMatrix = matrix(1, 1, 1, dimnames = list("a", "a")))
  expect_true(is.reversible(one_state))
})

test_that("is.reversible validates tolerance and rejects reducible chains", {
  expect_error(is.reversible(mc2, tolerance = -1), "tolerance")
  expect_error(is.reversible(mc2, tolerance = c(1e-8, 1e-8)), "tolerance")

  reducible <- new("markovchain", states = states2,
                   transitionMatrix = matrix(c(1, 0, 0.5, 0.5),
                                              byrow = TRUE, nrow = 2,
                                              dimnames = list(states2, states2)))
  expect_error(is.reversible(reducible), "irreducible")
})

## ---- mixingTime ------------------------------------------------------------

test_that("mixingTime finds the first t with TV distance below epsilon", {
  Pt <- P2
  d <- function(Pt) max(0.5 * rowSums(abs(sweep(Pt, 2, as.numeric(steadyStates(mc2)), "-"))))
  t <- 1L
  while (d(Pt) > 0.25) {
    Pt <- Pt %*% P2
    t <- t + 1L
  }
  expect_equal(mixingTime(mc2), t)
  expect_true(d(Pt) <= 0.25)
})

test_that("mixingTime is monotone non-decreasing as epsilon shrinks", {
  expect_lte(mixingTime(mc2, epsilon = 0.25), mixingTime(mc2, epsilon = 0.1))
  expect_lte(mixingTime(mc2, epsilon = 0.1), mixingTime(mc2, epsilon = 0.01))
})

test_that("mixingTime is zero for the one-state chain", {
  one_state <- new("markovchain", states = "a",
                   transitionMatrix = matrix(1, 1, 1, dimnames = list("a", "a")))
  expect_equal(mixingTime(one_state), 0L)
})

test_that("mixingTime supports column-stochastic storage", {
  mc2_col <- new("markovchain", states = states2, byrow = FALSE,
                 transitionMatrix = t(P2))
  expect_equal(mixingTime(mc2_col), mixingTime(mc2))
})

test_that("mixingTime rejects periodic and reducible chains, and validates arguments", {
  periodic <- new("markovchain", states = states2,
                  transitionMatrix = matrix(c(0, 1, 1, 0), 2, byrow = TRUE,
                                            dimnames = list(states2, states2)))
  expect_error(mixingTime(periodic), "aperiodic")

  reducible <- new("markovchain", states = states2,
                   transitionMatrix = matrix(c(1, 0, 0.5, 0.5),
                                              byrow = TRUE, nrow = 2,
                                              dimnames = list(states2, states2)))
  expect_error(mixingTime(reducible), "irreducible")

  expect_error(mixingTime(mc2, epsilon = 0), "epsilon")
  expect_error(mixingTime(mc2, epsilon = 1), "epsilon")
  expect_error(mixingTime(mc2, epsilon = -0.1), "epsilon")
  expect_error(mixingTime(mc2, maxIter = 0), "maxIter")
  expect_error(mixingTime(mc2, maxIter = 1.5), "maxIter")
})

test_that("mixingTime reports a clear error when maxIter is too small", {
  expect_error(mixingTime(mc2, epsilon = 1e-6, maxIter = 1L), "maxIter")
})

test_that("laziness restores a finite mixingTime to a periodic chain", {
  periodic <- new("markovchain", states = states2,
                  transitionMatrix = matrix(c(0, 1, 1, 0), 2, byrow = TRUE,
                                            dimnames = list(states2, states2)))
  lazyPeriodic <- lazyChain(periodic, alpha = 0.5)
  expect_equal(period(lazyPeriodic), 1L)
  expect_equal(mixingTime(lazyPeriodic), 1L) # already exactly stationary after 1 step
})
