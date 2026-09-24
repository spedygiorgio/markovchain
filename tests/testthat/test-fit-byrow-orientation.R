library(testthat)
library(markovchain)

# Regression tests for the byrow orientation of markovchainFit() results.
#
# Before this was fixed, every fitting path returned an *invalid* markovchain
# object when byrow = FALSE: mle/laplace/map transposed the transition matrix
# but left the byrow slot at its TRUE default, while bootstrap did the mirror
# image (slot set to FALSE, matrix left row-stochastic). Downstream methods
# that read the slot -- steadyStates(), is.irreducible(), slem(), ... -- then
# silently interpreted the matrix the wrong way round.

seqData <- c("a", "b", "c", "c", "e", "a", "b", "d", "e",
             "a", "c", "b", "e", "d", "a")
listData <- list(c("a", "b", "c", "c"), c("a", "c", "b", "c", "d"),
                 c("b", "d", "a", "c"))

test_that("every sequence fitting method returns a valid object for both orientations", {
  for (method in c("mle", "bootstrap", "laplace", "map")) {
    for (byrow in c(TRUE, FALSE)) {
      fit <- markovchainFit(seqData, method = method, byrow = byrow, nboot = 3)
      estimate <- fit$estimate

      expect_true(validObject(estimate),
                  info = paste(method, "byrow =", byrow))
      expect_equal(estimate@byrow, byrow,
                   info = paste(method, "byrow =", byrow))

      M <- estimate@transitionMatrix
      sums <- if (byrow) rowSums(M) else colSums(M)
      expect_equal(unname(sums), rep(1, nrow(M)), tolerance = 1e-10,
                   info = paste(method, "byrow =", byrow))
    }
  }
})

test_that("list input returns a valid object for both orientations", {
  for (method in c("mle", "map")) {
    for (byrow in c(TRUE, FALSE)) {
      estimate <- markovchainFit(listData, method = method, byrow = byrow)$estimate
      expect_true(validObject(estimate), info = paste(method, byrow))
      expect_equal(estimate@byrow, byrow, info = paste(method, byrow))
    }
  }
})

test_that("a byrow = FALSE fit is exactly the transpose of the byrow = TRUE fit", {
  for (method in c("mle", "laplace", "map")) {
    byRowFit <- markovchainFit(seqData, method = method, byrow = TRUE)
    byColFit <- markovchainFit(seqData, method = method, byrow = FALSE)
    expect_equal(byColFit$estimate@transitionMatrix,
                 t(byRowFit$estimate@transitionMatrix),
                 tolerance = 1e-12, info = method)
  }
})

test_that("uncertainty matrices follow the estimate's orientation", {
  # If the estimate is transposed but its standard errors are not, entry
  # (i, j) of the two would describe different transitions.
  for (method in c("mle", "map")) {
    byRowFit <- markovchainFit(seqData, method = method, byrow = TRUE)
    byColFit <- markovchainFit(seqData, method = method, byrow = FALSE)
    expect_equal(byColFit$standardError, t(byRowFit$standardError),
                 tolerance = 1e-12, info = method)
  }
})

test_that("downstream methods read both orientations identically", {
  byRowFit <- markovchainFit(seqData, method = "mle", byrow = TRUE)$estimate
  byColFit <- markovchainFit(seqData, method = "mle", byrow = FALSE)$estimate

  # steadyStates returns a row vector for row-stochastic storage and a column
  # vector for column-stochastic storage, so compare the values themselves.
  expect_equal(as.numeric(steadyStates(byRowFit)),
               as.numeric(steadyStates(byColFit)), tolerance = 1e-10)
  expect_equal(is.irreducible(byRowFit), is.irreducible(byColFit))
})

test_that("matrix input keeps byrow as a data-layout flag and returns row-stochastic estimates", {
  # For matrix/data.frame input, byrow says whether each trajectory is a row
  # or a column of the input data; the fitted chain is row-stochastic either
  # way. This is deliberately different from the sequence/list case.
  trajectories <- matrix(c("a", "b", "c",
                           "a", "c", "b"), nrow = 2, byrow = TRUE)

  for (byrow in c(TRUE, FALSE)) {
    fit <- markovchainFit(trajectories, method = "mle", byrow = byrow)
    expect_true(validObject(fit$estimate))
    expect_true(fit$estimate@byrow)
    expect_equal(unname(rowSums(fit$estimate@transitionMatrix)),
                 rep(1, nrow(fit$estimate@transitionMatrix)), tolerance = 1e-10)
    # The uncertainty matrices stay aligned with that row-stochastic estimate.
    expect_equal(dimnames(fit$standardError),
                 dimnames(fit$estimate@transitionMatrix))
  }
})
