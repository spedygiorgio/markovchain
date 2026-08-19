library(testthat)
library(markovchain)

states <- c("a", "b")
mc <- new("markovchain",
          states = states,
          transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
                                     byrow = TRUE, nrow = 2,
                                     dimnames = list(states, states)))

test_that("entropy rate matches the stationary definition", {
  pi <- c(0.25, 0.75)
  P <- mc@transitionMatrix
  logP <- matrix(0, nrow(P), ncol(P))
  positive <- P > 0
  logP[positive] <- log(P[positive], base = 2)
  expected <- -sum(pi * rowSums(P * logP))
  expect_equal(entropyRate(mc), expected, tolerance = 1e-12)
})

test_that("entropy rate supports alternative logarithm bases", {
  expect_equal(entropyRate(mc, base = exp(1)),
               entropyRate(mc, base = 2) * log(2),
               tolerance = 1e-12)
})

test_that("zero transition probabilities contribute zero", {
  # This irreducible chain contains a structural zero.
  structural_zero <- new("markovchain",
                         states = states,
                         transitionMatrix = matrix(c(0.5, 0.5, 1, 0),
                                                    byrow = TRUE, nrow = 2,
                                                    dimnames = list(states, states)))
  expect_true(is.finite(entropyRate(structural_zero)))
  expect_equal(entropyRate(structural_zero), 2 / 3, tolerance = 1e-12)
})

test_that("entropy rate is zero for a one-state chain", {
  one_state <- new("markovchain",
                   states = "a",
                   transitionMatrix = matrix(1, 1, 1,
                                             dimnames = list("a", "a")))
  expect_equal(entropyRate(one_state), 0, tolerance = 1e-12)
})

test_that("entropy rate validates the logarithm base", {
  expect_error(entropyRate(mc, base = 0), "base must be")
  expect_error(entropyRate(mc, base = 1), "base must be")
  expect_error(entropyRate(mc, base = c(2, 10)), "base must be")
})
