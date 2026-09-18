context("ctmcFit - lambda confidence interval")

### Regression test for the fix in ctmcFit() (src/ctmcFittingFunctions.cpp):
### the confidence interval for each state's exponential holding-time rate
### lambda_i had two bugs, both found by direct calculation:
###   1. the z quantile used was the ONE-sided normal quantile (1.645 for
###      confidencelevel = 0.95) instead of the 1.96 a two-sided interval
###      at that confidence level requires;
###   2. the lower and upper bounds were both derived from the very same
###      "factor" value (lambda_hat * (1 - z/sqrt(n))), just clipped
###      differently (max(0, factor) vs min(1, factor)) -- whenever factor
###      fell in [0, 1], the common case, this made lower and upper
###      IDENTICAL, i.e. a zero-width, non-informative interval. lambda is
###      a rate, not a probability, so it also should never be clipped
###      above at 1.

test_that("lambda confidence interval is not degenerate (lower != upper)", {
  data <- list(
    c("a", "b", "c", "a", "b", "a", "c", "b", "c"),
    c(0, 0.8, 2.1, 2.4, 4, 5, 5.9, 8.2, 9)
  )
  fit <- ctmcFit(data, confidencelevel = 0.95)

  lower <- fit$errors$lambdaConfidenceInterval$lowerEndpointVector
  upper <- fit$errors$lambdaConfidenceInterval$upperEndpointVector

  expect_true(all(upper > lower))
})

test_that("lambda confidence interval uses the correct two-sided z quantile", {
  data <- list(
    c("a", "b", "c", "a", "b", "a", "c", "b", "c"),
    c(0, 0.8, 2.1, 2.4, 4, 5, 5.9, 8.2, 9)
  )
  fit <- ctmcFit(data, confidencelevel = 0.95)

  stateData <- data[[1]]
  transData <- data[[2]]
  sortedStates <- sort(unique(stateData))
  stateCount <- setNames(numeric(length(sortedStates)), sortedStates)
  stateSojournTime <- setNames(numeric(length(sortedStates)), sortedStates)
  for (i in seq_len(length(stateData) - 1)) {
    stateCount[stateData[i]] <- stateCount[stateData[i]] + 1
    stateSojournTime[stateData[i]] <- stateSojournTime[stateData[i]] +
      (transData[i + 1] - transData[i])
  }

  z <- qnorm(1 - (1 - 0.95) / 2)  # the correct two-sided quantile, 1.96
  lambdaHat <- stateCount / stateSojournTime
  margin <- z * lambdaHat / sqrt(stateCount)
  expectedLower <- pmax(0, lambdaHat - margin)
  expectedUpper <- lambdaHat + margin

  lower <- fit$errors$lambdaConfidenceInterval$lowerEndpointVector
  upper <- fit$errors$lambdaConfidenceInterval$upperEndpointVector

  expect_equal(unname(lower), unname(expectedLower), tolerance = 1e-6)
  expect_equal(unname(upper), unname(expectedUpper), tolerance = 1e-6)
})

test_that("lambda upper confidence bound is not clipped at 1 (lambda is a rate, not a probability)", {
  # A short, fast-transitioning chain pushes lambda_hat well above 1.
  data <- list(
    c("a", "b", "a", "b", "a", "b", "a", "b", "a", "b"),
    seq(0, by = 0.05, length.out = 10)
  )
  fit <- ctmcFit(data, confidencelevel = 0.95)
  upper <- fit$errors$lambdaConfidenceInterval$upperEndpointVector

  expect_true(any(upper > 1))
})
