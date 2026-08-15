test_that("verifyEmpiricalToTheoretical has correct multinomial null", {
  P <- matrix(c(.5, .5, .5, .5), 2, 2, byrow = TRUE,
              dimnames = list(c("a", "b"), c("a", "b")))
  mc <- new("markovchain", states = c("a", "b"), transitionMatrix = P)
  counts <- matrix(c(50, 50, 50, 50), 2, 2, byrow = TRUE,
                   dimnames = dimnames(P))

  ans <- verifyEmpiricalToTheoretical(counts, mc, verbose = FALSE)
  expect_equal(ans$statistic, 0)
  expect_equal(ans$dof, 2)
  expect_equal(ans$p.value, 1)

  ans2 <- verifyEmpiricalToTheoretical(counts, mc, method = "Pearson", verbose = FALSE)
  expect_equal(ans2$statistic, 0)
})

test_that("verifyEmpiricalToTheoretical detects structural zeros", {
  P <- matrix(c(1, 0, 0, 1), 2, 2, byrow = TRUE,
              dimnames = list(c("a", "b"), c("a", "b")))
  mc <- new("markovchain", states = c("a", "b"), transitionMatrix = P)
  counts <- matrix(c(9, 1, 0, 10), 2, 2, byrow = TRUE,
                   dimnames = dimnames(P))
  ans <- verifyEmpiricalToTheoretical(counts, mc, verbose = FALSE)
  expect_true(is.infinite(ans$statistic))
  expect_equal(ans$p.value, 0)
})

test_that("verifyHomogeneity accepts sequences and matrices", {
  P <- matrix(c(.7, .3, .2, .8), 2, 2, byrow = TRUE,
              dimnames = list(c("a", "b"), c("a", "b")))
  mc <- new("markovchain", states = c("a", "b"), transitionMatrix = P)
  set.seed(123)
  x <- rmarkovchain(500, mc, t0 = "a")
  y <- rmarkovchain(500, mc, t0 = "b")
  ans <- verifyHomogeneity(list(x, y), verbose = FALSE)
  expect_true(is.finite(ans$statistic))
  expect_true(ans$dof > 0)
  expect_true(ans$p.value >= 0 && ans$p.value <= 1)
})

test_that("verifyMarkovProperty returns the revised result structure", {
  set.seed(123)
  P <- matrix(c(.7, .3, .2, .8), 2, 2, byrow = TRUE,
              dimnames = list(c("a", "b"), c("a", "b")))
  mc <- new("markovchain", states = c("a", "b"), transitionMatrix = P)
  x <- rmarkovchain(200, mc, t0 = "a")
  ans <- verifyMarkovProperty(x, method = "G", verbose = FALSE)
  expect_true(all(c("statistic", "dof", "p.value", "observed", "expected") %in% names(ans)))
  expect_true(is.finite(ans$statistic))
})
