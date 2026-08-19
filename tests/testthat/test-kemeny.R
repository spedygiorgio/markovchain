library(testthat)
library(markovchain)

states <- c("a", "b")
mc <- new("markovchain",
          states = states,
          transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
                                     byrow = TRUE, nrow = 2,
                                     dimnames = list(states, states)))

test_that("kemeny's constant matches the first-passage definition", {
  # pi = (0.25, 0.75), m_ab = 1 / 0.3 and m_ba = 1 / 0.1.
  expected <- 0.25 * (1 / 0.3) + 0.75 * (1 / 0.1)
  expect_equal(kemenyConstant(mc), expected, tolerance = 1e-12)
})

test_that("kemeny's constant is invariant to the starting state", {
  P <- mc@transitionMatrix
  pi <- as.numeric(steadyStates(mc))
  m <- matrix(0, 2, 2)
  m[1, 2] <- 1 / P[1, 2]
  m[2, 1] <- 1 / P[2, 1]
  m[1, 1] <- 1 / pi[1]
  m[2, 2] <- 1 / pi[2]
  values <- as.numeric(m %*% pi)
  expect_equal(values[1], values[2], tolerance = 1e-12)
  expect_equal(kemenyConstant(mc), values[1], tolerance = 1e-12)
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
