library(testthat)
library(markovchain)

states <- c("a", "b")
P <- matrix(c(0.7, 0.3, 0.1, 0.9),
            byrow = TRUE, nrow = 2,
            dimnames = list(states, states))
mc <- new("markovchain", states = states, transitionMatrix = P)

test_that("kemeny's constant matches the zero-diagonal hitting-time definition", {
  # pi = (0.25, 0.75), m_ab = 1/0.3 and m_ba = 1/0.1.
  # Starting from a: K = pi_b*m_ab; starting from b: K = pi_a*m_ba.
  expected_from_a <- 0.75 * (1 / 0.3)
  expected_from_b <- 0.25 * (1 / 0.1)
  expect_equal(expected_from_a, expected_from_b, tolerance = 1e-12)
  expect_equal(kemenyConstant(mc), expected_from_a, tolerance = 1e-12)
})

test_that("kemeny's constant is invariant to the starting state", {
  pi <- as.numeric(steadyStates(mc))
  m <- meanFirstPassageTime(mc)
  values <- as.numeric(m %*% pi)

  expect_equal(unname(diag(m)), c(0, 0), tolerance = 1e-12)
  expect_equal(values, rep(values[1], length(values)), tolerance = 1e-12)
  expect_equal(kemenyConstant(mc), values[1], tolerance = 1e-12)
})

test_that("the first-return convention differs by one", {
  pi <- as.numeric(steadyStates(mc))
  m_return <- meanFirstPassageTime(mc)
  diag(m_return) <- 1 / pi
  return_values <- as.numeric(m_return %*% pi)

  expect_equal(return_values, rep(kemenyConstant(mc) + 1, 2),
               tolerance = 1e-12)
})

test_that("kemeny's constant supports column-stochastic storage", {
  mc_col <- new("markovchain",
                states = states,
                byrow = FALSE,
                transitionMatrix = t(P))
  expect_equal(kemenyConstant(mc_col), kemenyConstant(mc), tolerance = 1e-12)
})

test_that("kemeny's constant is defined for periodic irreducible chains", {
  periodic <- new("markovchain",
                  states = states,
                  transitionMatrix = matrix(c(0, 1, 1, 0), 2, byrow = TRUE,
                                            dimnames = list(states, states)))
  expect_equal(kemenyConstant(periodic), 0.5, tolerance = 1e-12)
})

test_that("kemeny's constant agrees with the non-unit eigenvalue formula", {
  lambda <- eigen(P, only.values = TRUE)$values
  unit <- which.min(Mod(lambda - 1))
  spectral <- Re(sum(1 / (1 - lambda[-unit])))
  expect_equal(kemenyConstant(mc), spectral, tolerance = 1e-12)
})

test_that("kemeny's constant is zero for the one-state chain", {
  one_state <- new("markovchain",
                   states = "a",
                   transitionMatrix = matrix(1, 1, 1,
                                             dimnames = list("a", "a")))
  expect_equal(kemenyConstant(one_state), 0, tolerance = 1e-12)
})

test_that("kemeny's constant rejects reducible chains", {
  reducible <- new("markovchain",
                   states = states,
                   transitionMatrix = matrix(c(1, 0, 0.5, 0.5),
                                              byrow = TRUE, nrow = 2,
                                              dimnames = list(states, states)))
  expect_error(kemenyConstant(reducible), "defined here only for irreducible")
})
