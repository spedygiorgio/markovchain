context("predictHommc")

# two sequences on states a, b; matrices stored by column (P[to, from]) as in
# fitHighOrderMultivarMC(). Index of P_h^{(jk)}: n * s * (j - 1) + n * (k - 1) + h
.twoSequenceModel <- function(lambda, P, byrow = FALSE) {
  new("hommc", order = 1, states = c("a", "b"), P = P, Lambda = lambda,
      byrow = byrow, name = "test")
}

test_that("predictHommc uses the past of the source sequence and reads P by column", {
  P <- array(0, dim = c(2, 2, 4), dimnames = list(c("a", "b"), c("a", "b")))
  # sequence 1 is driven only by sequence 2 (P^{(12)} = index 2): from a -> b, from b -> b
  P[, , 2] <- matrix(c(0, 1,
                       0, 1), 2)          # columns: from a, from b
  P[, , 1] <- diag(2)
  # sequence 2 is driven only by itself (P^{(22)} = index 4): from a -> a, from b -> a
  P[, , 4] <- matrix(c(1, 0,
                       1, 0), 2)
  P[, , 3] <- diag(2)
  model <- .twoSequenceModel(c(0, 1, 0, 1), P)
  init <- matrix(c("a", "a"), nrow = 2)
  out <- predictHommc(model, 3, init)
  # sequence 1 at each step: P^{(12)} applied to the previous state of sequence 2 -> always b
  expect_equal(out[1, ], c("b", "b", "b"))
  expect_equal(out[2, ], c("a", "a", "a"))

  # the same model stored by row gives the same result
  Prow <- P
  for (i in 1:4) Prow[, , i] <- t(P[, , i])
  expect_equal(predictHommc(.twoSequenceModel(c(0, 1, 0, 1), Prow, byrow = TRUE), 3, init), out)
})

test_that("all sequences are drawn from the same past", {
  P <- array(0, dim = c(2, 2, 4), dimnames = list(c("a", "b"), c("a", "b")))
  swap <- matrix(c(0, 1, 1, 0), 2)
  P[, , 1] <- swap          # sequence 1 from itself
  P[, , 2] <- diag(2)
  P[, , 3] <- diag(2)       # sequence 2 copies sequence 1
  P[, , 4] <- diag(2)
  model <- .twoSequenceModel(c(1, 0, 1, 0), P)
  out <- predictHommc(model, 4, matrix(c("a", "b"), nrow = 2))
  expect_equal(out[1, ], c("b", "a", "b", "a"))
  # sequence 2 copies the previous state of sequence 1, not the one just drawn
  expect_equal(out[2, ], c("a", "b", "a", "b"))
})

test_that("predictHommc works with order 2 and reproduces the model probabilities", {
  set.seed(1)
  P <- array(0, dim = c(2, 2, 2), dimnames = list(c("a", "b"), c("a", "b")))
  P[, , 1] <- matrix(c(0.9, 0.1, 0.2, 0.8), 2)  # lag 1, by column
  P[, , 2] <- matrix(c(0.5, 0.5, 0.5, 0.5), 2)  # lag 2
  model <- new("hommc", order = 2, states = c("a", "b"), P = P,
               Lambda = c(0.6, 0.4), byrow = FALSE, name = "univariate")
  draws <- replicate(4000, predictHommc(model, 1, matrix(c("b", "a"), nrow = 1))[1, 1])
  # P(next = a | last = a, previous = b) = 0.6 * 0.9 + 0.4 * 0.5 = 0.74
  expect_equal(mean(draws == "a"), 0.74, tolerance = 0.03)
  expect_length(predictHommc(model, 5, matrix(c("b", "a"), nrow = 1)), 5)
})

test_that("predictHommc validates its arguments", {
  P <- array(diag(2), dim = c(2, 2, 4), dimnames = list(c("a", "b"), c("a", "b")))
  model <- .twoSequenceModel(rep(0.5, 4), P)
  expect_error(predictHommc(model, 0, matrix("a", 2, 1)), "positive integer")
  expect_error(predictHommc(model, 2, matrix("a", 1, 1)), "previous states")
  expect_error(predictHommc(model, 2, matrix(c("a", "z"), 2, 1)), "invalid states")
  expect_error(predictHommc(list(), 2), "hommc")
})
