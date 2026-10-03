context("markovchainFit progress bar (#232)")

seqProgress <- c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a",
                 "b", "b", "b", "a", "c", "b", "c", "a")

test_that("the bootstrap shows a progress bar without changing the results", {
  set.seed(10)
  quiet <- markovchainFit(seqProgress, method = "bootstrap", nboot = 20)
  set.seed(10)
  out <- capture.output(shown <- markovchainFit(seqProgress, method = "bootstrap",
                                                nboot = 20, progress = TRUE))
  expect_true(any(grepl("100%", out)))
  expect_true(any(grepl("=", out)))
  expect_equal(shown$estimate, quiet$estimate)
  expect_equal(shown$standardError, quiet$standardError)
})

test_that("the progress bar also works with the parallel bootstrap", {
  out <- capture.output(fit <- markovchainFit(seqProgress, method = "bootstrap", nboot = 10,
                                              parallel = TRUE, progress = TRUE))
  expect_true(any(grepl("100%", out)))
  expect_s4_class(fit$estimate, "markovchain")
})

test_that("progress is ignored by the fast methods and validated", {
  for (method in c("mle", "laplace", "map")) {
    out <- capture.output(fit <- markovchainFit(seqProgress, method = method, progress = TRUE))
    expect_length(out, 0)
    expect_equal(fit, markovchainFit(seqProgress, method = method))
  }
  expect_error(markovchainFit(seqProgress, progress = NA), "progress")
  expect_error(markovchainFit(seqProgress, progress = "yes"), "progress")
})
