context("C++ entry points reject inputs that used to read or write out of bounds")

# Each of these used to crash R (segfault), read past the end of a vector, or
# loop forever. They must now fail with an ordinary R error, or return FALSE.

test_that("ctmcFit checks the shape of its data", {
  expect_error(ctmcFit(list(c("a", "b", "c", "a"))), "list with the visited states")
  expect_error(ctmcFit(list(c("a", "b", "c", "a"), c(0, 1))), "same length")
  ok <- ctmcFit(list(c("a", "b", "a", "b"), c(0, 1, 2, 3)))
  expect_type(ok, "list")
})

test_that("generatorToTransitionMatrix requires a square generator", {
  expect_error(generatorToTransitionMatrix(matrix(c(-1, 1, 2, -2, 0, 0), 3, 2, byrow = TRUE)),
               "square")
  expect_error(generatorToTransitionMatrix(matrix(c(-1, 1, 0, 2, -2, 0), 2, 3, byrow = TRUE),
                                           byrow = FALSE), "square")
  gen <- matrix(c(-1, 1, 2, -2), 2, byrow = TRUE)
  expect_equal(unname(rowSums(generatorToTransitionMatrix(gen))), c(1, 1))
})

test_that("multinomial confidence intervals need counts shaped like the matrix", {
  expect_error(markovchain:::multinomialConfidenceIntervals(
    matrix(.5, 2, 2), matrix(c(3, 4, 5, 1, 2, 3), 2, 3), .95), "same dimensions")
})

test_that("MAP fit validates the hyperparameter matrix", {
  s <- c("a", "b", "a", "b", "b")
  expect_error(markovchainFit(s, method = "map", hyperparam = matrix(1, 2, 2)),
               "row names and column names")
  hp <- matrix(c(1, 1, 1, NA), 2, dimnames = list(c("a", "b"), c("a", "b")))
  expect_error(markovchainFit(s, method = "map", hyperparam = hp), "greater than or equal to 1")
  hp[2, 2] <- NaN
  expect_error(markovchainFit(s, method = "map", hyperparam = hp), "greater than or equal to 1")
})

test_that("NaN in a transition matrix is rejected instead of looping", {
  expect_false(markovchain:::.isStochasticMatrix(matrix(c(NaN, 1, .5, .5), 2, byrow = TRUE), TRUE))
  expect_false(markovchain:::.isStochasticMatrix(matrix(c(NA, NA, .5, .5), 2, byrow = TRUE), TRUE))
  expect_error(new("markovchain",
                   transitionMatrix = matrix(c(NA, NA, .5, .5), 2, byrow = TRUE,
                                             dimnames = list(c("a", "b"), c("a", "b")))))
})

test_that("the parallel simulation refuses a markovchain whose slots disagree", {
  P <- matrix(.5, 2, 2, dimnames = list(c("a", "b"), c("a", "b")))
  mc <- new("markovchain", states = c("a", "b"), transitionMatrix = P)
  mc@states <- "a"   # slot assignment bypasses validity
  ml <- new("markovchainList", markovchains = list(mc))
  expect_error(markovchain:::.markovchainSequenceParallelRcpp(ml, 5L, FALSE, character()),
               "consistent dimensions")
})

test_that("the imprecise-probability kernel checks its indices", {
  states <- c("n", "y")
  Q <- matrix(c(-1, 1, 1, -1), 2, byrow = TRUE, dimnames = list(states, states))
  range <- matrix(c(1/52, 3/52, 1/2, 2), 2, byrow = TRUE)
  ic <- new("ictmc", states = states, Q = Q, range = range, name = "x")
  expect_error(markovchain:::.impreciseProbabilityatTRCpp(ic, 0L, 0L, 1L, 0.1), "between 1 and")
  expect_error(markovchain:::.impreciseProbabilityatTRCpp(ic, 3L, 0L, 1L, 0.1), "between 1 and")
  ok <- markovchain:::.impreciseProbabilityatTRCpp(ic, 1L, 0L, 1L, 0.1)
  expect_length(ok, 2L)
})

test_that("a matrix with negative entries is not accepted as stochastic", {
  # the row sums to 1, but 1.5 and -0.5 are not probabilities
  expect_false(markovchain:::.isStochasticMatrix(matrix(c(1.5, -.5, .5, .5), 2, byrow = TRUE), TRUE))
  expect_false(markovchain:::.isStochasticMatrix(matrix(c(1.5, .5, -.5, .5), 2, byrow = FALSE), FALSE))
  expect_true(markovchain:::.isStochasticMatrix(matrix(c(.5, .5, .4, .6), 2, byrow = TRUE), TRUE))
})

test_that("ctmcFit(byrow = FALSE) returns a column-oriented generator", {
  d <- list(c("a", "b", "a", "c", "a", "b", "c", "a"), c(0, 1, 2.5, 3, 4.5, 5, 7, 8))
  byRow <- ctmcFit(d, byrow = TRUE)$estimate
  byCol <- ctmcFit(d, byrow = FALSE)$estimate
  expect_true(byRow@byrow)
  expect_false(byCol@byrow)
  # rows of one, columns of the other sum to zero: a proper generator
  expect_equal(unname(rowSums(byRow@generator)), rep(0, 3), tolerance = 1e-12)
  expect_equal(unname(colSums(byCol@generator)), rep(0, 3), tolerance = 1e-12)
  # and the two are transposes of each other
  expect_equal(unname(byCol@generator), unname(t(byRow@generator)), tolerance = 1e-12)
})
