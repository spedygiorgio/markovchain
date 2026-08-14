test_that("is.stochasticallyMonotone recognizes monotone transition matrices", {
  monotoneMatrix <- matrix(c(0.8, 0.2,
                             0.3, 0.7), nrow = 2, byrow = TRUE)
  monotoneChain <- as(monotoneMatrix, "markovchain")

  expect_true(is.stochasticallyMonotone(monotoneMatrix))
  expect_true(is.stochasticallyMonotone(monotoneChain))
  expect_true(is.stochasticallyMonotone(matrix(1, nrow = 1)))
})

test_that("is.stochasticallyMonotone rejects non-monotone matrices and invalid inputs", {
  nonMonotoneMatrix <- matrix(c(0.3, 0.7,
                                0.8, 0.2), nrow = 2, byrow = TRUE)

  expect_false(is.stochasticallyMonotone(nonMonotoneMatrix))
  expect_error(is.stochasticallyMonotone(1:2), "must be a `markovchain` object or a matrix")
})