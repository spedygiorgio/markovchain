context("higher order fits: input types, small samples, missing Rsolnp")

.seqHO <- function() {
  c("a", "a", "b", "b", "a", "c", "b", "a", "b", "c", "a", "b",
    "c", "a", "b", "c", "a", "b", "a", "b")
}

test_that("fitHigherOrder accepts factors and numbers like fitMTD and higherOrderLogLik", {
  s <- .seqHO()
  ref <- fitHigherOrder(s, 2, method = "mle")
  expect_equal(fitHigherOrder(factor(s), 2, method = "mle"), ref)
  num <- c(1, 2, 3, 1, 2, 3, 3, 2, 1, 1, 2, 3, 2, 1, 3, 1, 2, 2, 3, 1)
  expect_equal(fitHigherOrder(num, 2, method = "mle"),
               fitHigherOrder(as.character(num), 2, method = "mle"))
  expect_equal(rownames(fitHigherOrder(num, 1, method = "mle")$Q[[1L]]), c("1", "2", "3"))
})

test_that("fitHigherOrder still rejects unusable sequences", {
  expect_error(fitHigherOrder(c("a", NA, "b", "a"), 1, method = "mle"), "missing values")
  expect_error(fitHigherOrder("a", 1, method = "mle"), "at least two")
  expect_error(fitHigherOrder(list(.seqHO()), 1, method = "mle"), "vector of at least two states")
  expect_error(fitHigherOrder(NULL, 1, method = "mle"), "vector of at least two states")
})

test_that("a warning is given when observations do not exceed the free parameters", {
  # 4 states in 4 observations, order 3: one observation, 4 free parameters
  expect_warning(fit <- fitHigherOrder(c("a", "b", "a", "b"), 3, method = "mle"),
                 "only 1 observation can be used to estimate the 4 free parameters")
  expect_equal(unname(fit$lambda), rep(1 / 3, 3))
  expect_warning(fitMTD(c("a", "b", "a", "b", "a"), 3), "only 2 observations can be used")
  # nothing to report on a sensible sample
  expect_no_warning(fitHigherOrder(.seqHO(), 2, method = "mle"))
  expect_no_warning(fitMTD(.seqHO(), 2))
})

test_that("the warning does not change the fit", {
  s <- c("a", "b", "a", "b", "a")
  withCallingHandlers(a <- fitMTD(s, 3), warning = function(w) invokeRestart("muffleWarning"))
  expect_equal(a$nobs, 2L)
  expect_equal(a$npar, 4L)
})

test_that("fitHigherOrder(method = 'lsq') warns and returns NULL when Rsolnp is unavailable", {
  testthat::local_mocked_bindings(requireNamespace = function(package, ..., quietly = FALSE) FALSE,
                                  .package = "base")
  expect_warning(out <- fitHigherOrder(.seqHO(), 2), "Rsolnp is required")
  expect_null(out)
  # the maximum likelihood fit does not need Rsolnp
  expect_no_warning(mle <- fitHigherOrder(.seqHO(), 2, method = "mle"))
  expect_length(mle$lambda, 2L)
  expect_warning(expect_error(higherOrderLogLik(.seqHO(), order = 2), "Rsolnp is required"),
                 "Rsolnp is required")
})
