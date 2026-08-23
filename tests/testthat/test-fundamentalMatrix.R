test_that("fundamentalMatrix returns the inverse complement of Q", {
  states <- c("a", "b", "absorbed")
  mc <- new("markovchain", states = states,
    transitionMatrix = matrix(c(
      0.5, 0.4, 0.1,
      0.2, 0.6, 0.2,
      0,   0,   1
    ), nrow = 3, byrow = TRUE,
    dimnames = list(states, states)))

  q <- matrix(c(0.5, 0.4, 0.2, 0.6), nrow = 2, byrow = TRUE)
  expected <- solve(diag(2) - q)

  expect_equal(fundamentalMatrix(mc), expected)
  expect_identical(rownames(fundamentalMatrix(mc)), c("a", "b"))
  expect_identical(colnames(fundamentalMatrix(mc)), c("a", "b"))
})

test_that("fundamentalMatrix rejects non-absorbing chains", {
  mc <- new("markovchain", states = c("a", "b"),
    transitionMatrix = matrix(c(0.7, 0.3, 0.2, 0.8),
      nrow = 2, byrow = TRUE,
      dimnames = list(c("a", "b"), c("a", "b"))))

  expect_error(fundamentalMatrix(mc), "absorbing Markov chain")
})

test_that("fundamentalMatrix handles a chain with no transient states", {
  states <- c("a", "b")
  mc <- new("markovchain", states = states,
    transitionMatrix = diag(2),
    dimnames = list(states, states))

  out <- fundamentalMatrix(mc)
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(0L, 0L))
})
