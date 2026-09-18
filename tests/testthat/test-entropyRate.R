library(testthat)
library(markovchain)

states <- c("a", "b")
P <- matrix(c(0.7, 0.3, 0.1, 0.9),
            byrow = TRUE, nrow = 2,
            dimnames = list(states, states))
mc <- new("markovchain", states = states, transitionMatrix = P)

test_that("entropy rate matches the stationary definition", {
  pi <- c(0.25, 0.75)
  expected <- -sum(pi * rowSums(P * log(P, base = 2)))
  expect_equal(entropyRate(mc), expected, tolerance = 1e-12)
})

test_that("entropy rate supports alternative logarithm bases", {
  expect_equal(entropyRate(mc, base = exp(1)),
               entropyRate(mc, base = 2) * log(2),
               tolerance = 1e-12)
  expect_equal(entropyRate(mc, base = 10),
               entropyRate(mc, base = 2) * log10(2),
               tolerance = 1e-12)
})

test_that("zero transition probabilities contribute zero", {
  structural_zero <- new("markovchain",
                         states = states,
                         transitionMatrix = matrix(c(0.5, 0.5, 1, 0),
                                                    byrow = TRUE, nrow = 2,
                                                    dimnames = list(states, states)))
  expect_true(is.finite(entropyRate(structural_zero)))
  expect_equal(entropyRate(structural_zero), 2 / 3, tolerance = 1e-12)
})

test_that("entropy rate supports column-stochastic storage", {
  mc_col <- new("markovchain", states = states, byrow = FALSE,
                transitionMatrix = t(P))
  expect_equal(entropyRate(mc_col), entropyRate(mc), tolerance = 1e-12)
})

test_that("a deterministic periodic chain has zero entropy rate", {
  periodic <- new("markovchain",
                  states = states,
                  transitionMatrix = matrix(c(0, 1, 1, 0), 2, byrow = TRUE,
                                            dimnames = list(states, states)))
  expect_equal(entropyRate(periodic), 0, tolerance = 1e-12)
})

test_that("an iid uniform chain has maximal entropy rate", {
  iid_states <- c("a", "b", "c")
  iid <- new("markovchain", states = iid_states,
             transitionMatrix = matrix(1 / 3, 3, 3,
                                       dimnames = list(iid_states, iid_states)))
  expect_equal(entropyRate(iid), log2(3), tolerance = 1e-12)
  expect_gte(entropyRate(iid), 0)
  expect_lte(entropyRate(iid), log2(length(iid_states)))
})

test_that("entropy rate is zero for a one-state chain", {
  one_state <- new("markovchain",
                   states = "a",
                   transitionMatrix = matrix(1, 1, 1,
                                             dimnames = list("a", "a")))
  expect_equal(entropyRate(one_state), 0, tolerance = 1e-12)
})

test_that("entropy rate validates the logarithm base", {
  for (invalid in list(-2, 0, 0.5, 1, Inf, NA_real_, c(2, 10), "2")) {
    expect_error(entropyRate(mc, base = invalid), "base must be")
  }
})

test_that("entropy rate rejects reducible chains", {
  reducible <- new("markovchain", states = states,
                   transitionMatrix = diag(2),
                   name = "reducible")
  expect_error(entropyRate(reducible), "irreducible")
})
