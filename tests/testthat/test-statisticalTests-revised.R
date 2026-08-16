test_that("verifyEmpiricalToTheoretical has correct multinomial null", {
  P <- matrix(c(.5, .5, .5, .5), 2, 2, byrow = TRUE,
              dimnames = list(c("a", "b"), c("a", "b")))
  mc <- new("markovchain", states = c("a", "b"), transitionMatrix = P)
  counts <- matrix(c(50, 50, 50, 50), 2, 2, byrow = TRUE,
                   dimnames = dimnames(P))

  ans <- verifyEmpiricalToTheoretical(counts, mc, verbose = FALSE)
  expect_s3_class(ans, "htest")
  expect_equal(unname(ans$statistic), 0)
  expect_equal(unname(ans$parameter), 2)
  expect_equal(ans$dof, 2)
  expect_equal(ans$p.value, 1)
  expect_equal(names(ans$statistic), "G-squared")
  expect_equal(names(ans$parameter), "df")

  ans2 <- verifyEmpiricalToTheoretical(counts, mc, method = "Pearson", verbose = FALSE)
  expect_s3_class(ans2, "htest")
  expect_equal(unname(ans2$statistic), 0)
  expect_equal(names(ans2$statistic), "X-squared")
})

test_that("verifyEmpiricalToTheoretical detects structural zeros", {
  P <- matrix(c(1, 0, 0, 1), 2, 2, byrow = TRUE,
              dimnames = list(c("a", "b"), c("a", "b")))
  mc <- new("markovchain", states = c("a", "b"), transitionMatrix = P)
  counts <- matrix(c(9, 1, 0, 10), 2, 2, byrow = TRUE,
                   dimnames = dimnames(P))
  ans <- verifyEmpiricalToTheoretical(counts, mc, verbose = FALSE)
  expect_s3_class(ans, "htest")
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
  expect_s3_class(ans, "htest")
  expect_true(is.finite(ans$statistic))
  expect_true(ans$dof > 0)
  expect_true(ans$parameter[["df"]] > 0)
  expect_true(ans$p.value >= 0 && ans$p.value <= 1)
})

test_that("verifyMarkovProperty returns the revised htest structure", {
  set.seed(123)
  P <- matrix(c(.7, .3, .2, .8), 2, 2, byrow = TRUE,
              dimnames = list(c("a", "b"), c("a", "b")))
  mc <- new("markovchain", states = c("a", "b"), transitionMatrix = P)
  x <- rmarkovchain(200, mc, t0 = "a")
  ans <- verifyMarkovProperty(x, method = "G", verbose = FALSE)
  expect_s3_class(ans, "htest")
  expect_true(all(c("statistic", "parameter", "p.value", "method", "data.name") %in% names(ans)))
  expect_true(all(c("observed", "expected", "transitionMatrix") %in% names(ans)))
  expect_true(is.finite(ans$statistic))
  expect_equal(names(ans$parameter), "df")
})

test_that("assessStationarity returns a standard htest object", {
  set.seed(123)
  P <- matrix(c(.7, .3, .2, .8), 2, 2, byrow = TRUE,
              dimnames = list(c("a", "b"), c("a", "b")))
  mc <- new("markovchain", states = c("a", "b"), transitionMatrix = P)
  x <- rmarkovchain(400, mc, t0 = "a")
  ans <- assessStationarity(x, nblocks = 4, verbose = FALSE)
  expect_s3_class(ans, "htest")
  expect_true(is.finite(ans$statistic))
  expect_true(ans$parameter[["df"]] > 0)
  expect_equal(names(ans$statistic), "X-squared")
  expect_true(grepl("time-homogeneity", ans$method, fixed = TRUE))

  printed <- capture.output(print(ans))
  expect_true(any(grepl("Pearson's Chi-squared test for time-homogeneity", printed, fixed = TRUE)))
  expect_true(any(grepl("X-squared", printed, fixed = TRUE)))
  expect_true(any(grepl("df", printed, fixed = TRUE)))
  expect_true(any(grepl("p-value", printed, fixed = TRUE)))
})
