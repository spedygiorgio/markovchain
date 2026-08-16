test_that("assessStationarity handles explicitly declared structural zeros", {
  # Alternating states ensures that 1 -> 1 is genuinely absent.
  sequence <- c(1, 2, 1, 2, 1, 2, 1, 2,
                1, 2, 1, 2, 1, 2, 1, 2)

  structural.zeros <- matrix(FALSE, 2, 2,
                              dimnames = list(c("1", "2"), c("1", "2")))
  structural.zeros[1, 1] <- TRUE

  ans <- assessStationarity(sequence, nblocks = 4,
                            structural.zeros = structural.zeros,
                            verbose = FALSE)

  expect_true(is.finite(ans$statistic))
  expect_true(is.finite(ans$dof))
  expect_true(is.finite(ans$p.value))
})

test_that("sampling zeros are retained in assessStationarity", {
  # 1 -> 2 is absent in the first block but is observed later; it is a
  # sampling zero, not a structural zero.
  sequence <- c(1, 1, 1, 1, 2, 1, 2, 2,
                1, 2, 1, 2, 1, 2, 1, 2)

  ans <- assessStationarity(sequence, nblocks = 2, verbose = FALSE)

  expect_true(is.finite(ans$statistic))
  expect_true(ans$dof > 0)
  expect_true(is.finite(ans$p.value))
})

test_that("structural zeros must not occur in the sequence", {
  sequence <- c(1, 2, 1, 2, 1, 2)
  structural.zeros <- matrix(FALSE, 2, 2)
  structural.zeros[1, 2] <- TRUE

  expect_error(
    assessStationarity(sequence, nblocks = 2,
                       structural.zeros = structural.zeros,
                       verbose = FALSE),
    "structural zero"
  )
})

test_that("named structural-zero matrices are matched to sequence states", {
  sequence <- c("b", "a", "b", "a", "b", "a", "b", "a")
  structural.zeros <- matrix(FALSE, 2, 2,
                             dimnames = list(c("a", "b"), c("a", "b")))
  structural.zeros["a", "a"] <- TRUE

  ans <- assessStationarity(sequence, nblocks = 2,
                            structural.zeros = structural.zeros,
                            verbose = FALSE)

  expect_true(is.finite(ans$statistic))
  expect_true(ans$dof > 0)
})

test_that("issue 63 regression: structural zeros no longer produce NaN", {
  # Reduced version of the pattern reported in issue #63: state 1 never
  # transitions to state 3, so the corresponding cell is structurally zero.
  sequence <- c(2, 5, 2, 4, 2, 5, 2, 4, 2, 4, 2, 2, 5, 2, 4, 2,
                5, 3, 4, 3, 4, 1, 4, 4, 1, 2, 5, 2, 1, 5, 4, 2,
                5, 1, 4, 4, 2, 3, 6, 6, 2, 2, 6, 5, 5, 5, 6, 3,
                6, 3, 3, 6, 3, 2, 5, 3, 5, 2, 6, 2, 2, 2, 5, 2,
                5, 2, 5, 5, 4, 2, 5, 5, 5, 2, 5, 2, 5, 2, 5, 2,
                2, 5, 2, 5, 2, 2, 5, 5, 2, 2, 2, 5, 2, 5, 2, 5,
                2, 5, 2, 2, 5, 2, 5, 3, 5, 2, 5, 5, 2, 2, 5, 2,
                2, 5, 2, 2, 5, 2, 6, 2, 5, 4, 2, 5, 2, 4, 4, 2,
                1, 2, 5, 2, 5, 2, 5, 5, 5, 3, 5, 2, 5, 5, 2, 5,
                5, 2, 5, 5, 2, 5, 2, 5, 5, 2, 1, 2, 5, 1, 5, 1,
                5, 5, 4, 2, 2, 2, 3, 6, 3, 6, 3, 6, 3, 6, 3, 6)

  structural.zeros <- matrix(FALSE, 6, 6,
                             dimnames = list(as.character(1:6), as.character(1:6)))
  structural.zeros["1", "3"] <- TRUE

  ans <- assessStationarity(sequence, nblocks = 8,
                            structural.zeros = structural.zeros,
                            verbose = FALSE)

  expect_true(is.finite(ans$statistic))
  expect_true(is.finite(ans$dof))
  expect_true(is.finite(ans$p.value))
})
