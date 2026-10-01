context("assessIndependence")

test_that("the Pearson statistic coincides with chisq.test on the transition count table", {
  set.seed(123)
  x <- sample(c("a", "b", "c"), 300, replace = TRUE)
  res <- assessIndependence(x, verbose = FALSE)
  ref <- suppressWarnings(chisq.test(createSequenceMatrix(x), correct = FALSE))

  expect_s3_class(res, "htest")
  expect_equal(unname(res$statistic), unname(ref$statistic))
  expect_equal(unname(res$parameter), unname(ref$parameter))
  expect_equal(res$p.value, ref$p.value)
})

test_that("the G statistic is the likelihood-ratio statistic of the count table", {
  set.seed(7)
  x <- sample(c("a", "b", "c"), 300, replace = TRUE)
  res <- assessIndependence(x, method = "G", verbose = FALSE)
  O <- res$observed
  E <- res$expected
  expect_equal(unname(res$statistic), 2 * sum(O[O > 0] * log(O[O > 0] / E[O > 0])))
  expect_equal(unname(res$parameter), (nrow(O) - 1) * (ncol(O) - 1))
})

test_that("counts, expected counts and degrees of freedom are consistent", {
  x <- c("a", "b", "a", "a", "b", "b", "a", "b", "a", "a", "b", "a")
  res <- assessIndependence(x, verbose = FALSE)
  expect_equal(sum(res$observed), length(x) - 1)
  expect_equal(sum(res$expected), length(x) - 1)
  expect_equal(unname(rowSums(res$expected)), unname(rowSums(res$observed)))
  expect_equal(unname(colSums(res$expected)), unname(colSums(res$observed)))
  expect_equal(unname(res$parameter), 1)
})

test_that("an independent sequence is not rejected and a Markov one is", {
  set.seed(2024)
  iid <- sample(c("a", "b", "c"), 2000, replace = TRUE)
  expect_gt(assessIndependence(iid, verbose = FALSE)$p.value, 0.01)

  mc <- new("markovchain", states = c("a", "b"),
            transitionMatrix = matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE))
  dep <- rmarkovchain(500, mc, t0 = "a")
  expect_lt(assessIndependence(dep, verbose = FALSE)$p.value, 1e-6)
  expect_lt(assessIndependence(dep, method = "G", verbose = FALSE)$p.value, 1e-6)
})

test_that("states never used as departure or arrival do not enter the degrees of freedom", {
  # "c" occurs only as the very last observation: an arrival state (a column)
  # but never a departure state (a row)
  x <- c("a", "b", "a", "b", "a", "b", "a", "b", "a", "c")
  res <- assessIndependence(x, verbose = FALSE)
  expect_equal(dim(res$observed), c(2L, 3L))
  expect_equal(unname(res$parameter), (2 - 1) * (3 - 1))
})

test_that("degenerate sequences give zero degrees of freedom and an NA p-value", {
  res <- assessIndependence(rep("a", 10), verbose = FALSE)
  expect_equal(unname(res$parameter), 0)
  expect_true(is.na(res$p.value))
})

test_that("factor and numeric sequences are accepted", {
  set.seed(1)
  x <- sample(1:3, 200, replace = TRUE)
  expect_equal(assessIndependence(x, verbose = FALSE)$statistic,
               assessIndependence(as.character(x), verbose = FALSE)$statistic)
  expect_equal(assessIndependence(factor(x), verbose = FALSE)$statistic,
               assessIndependence(as.character(x), verbose = FALSE)$statistic)
})

test_that("verbose controls printing and the result is returned invisibly", {
  x <- c("a", "b", "a", "a", "b", "b", "a", "b", "a", "a", "b", "a")
  expect_output(assessIndependence(x, verbose = TRUE), "independence")
  expect_silent(assessIndependence(x, verbose = FALSE))
  expect_invisible(assessIndependence(x, verbose = FALSE))
})

test_that("invalid input is rejected", {
  expect_error(assessIndependence(c("a", "b")), "at least three")
  expect_error(assessIndependence(c("a", "b", NA, "a")), "missing")
  expect_error(assessIndependence(c("a", "b", "a", "b"), method = "simulation"))
})
