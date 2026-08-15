test_that("verifyMarkovProperty returns the correct df for a binary chain", {
  sequence <- c(
    "a", "b", "a", "a", "a", "a", "b", "a", "b",
    "a", "b", "a", "a", "b", "b", "b", "a"
  )

  ans <- verifyMarkovProperty(
    sequence,
    method = "Pearson",
    verbose = FALSE
  )

  expect_equal(ans$dof, 2)
  expect_true(is.finite(ans$statistic))
  expect_true(is.finite(ans$p.value))
})


test_that("G and Pearson statistics are non-negative", {
  sequence <- c(
    "a", "b", "a", "a", "a", "a", "b", "a", "b",
    "a", "b", "a", "a", "b", "b", "b", "a"
  )

  g <- verifyMarkovProperty(sequence, method = "G", verbose = FALSE)
  p <- verifyMarkovProperty(sequence, method = "Pearson", verbose = FALSE)

  expect_gte(g$statistic, 0)
  expect_gte(p$statistic, 0)
  expect_equal(g$dof, p$dof)
})


test_that("simulation method is reproducible with a seed", {
  sequence <- c(
    "a", "b", "a", "a", "a", "a", "b", "a", "b",
    "a", "b", "a", "a", "b", "b", "b", "a"
  )

  ans1 <- verifyMarkovProperty(
    sequence,
    method = "simulation",
    B = 99,
    seed = 123,
    verbose = FALSE
  )
  ans2 <- verifyMarkovProperty(
    sequence,
    method = "simulation",
    B = 99,
    seed = 123,
    verbose = FALSE
  )

  expect_equal(ans1$p.value, ans2$p.value)
  expect_equal(ans1$simulations, ans2$simulations)
})


test_that("short sequences are rejected", {
  expect_error(
    verifyMarkovProperty(c("a", "b", "a"), verbose = FALSE),
    "at least four"
  )
})


test_that("missing values are rejected", {
  expect_error(
    verifyMarkovProperty(c("a", NA, "b", "a"), verbose = FALSE),
    "missing"
  )
})


test_that("simulation requires an outgoing transition for every state", {
  sequence <- c("a", "b", "b", "b", "c")

  expect_error(
    verifyMarkovProperty(
      sequence,
      method = "simulation",
      B = 9,
      verbose = FALSE
    ),
    "outgoing transition"
  )
})
