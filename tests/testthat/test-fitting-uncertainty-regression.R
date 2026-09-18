test_that("MLE confidence intervals do not degenerate for zero or full observed cells", {
  seq <- c(rep(c("a", "a"), 20), rep(c("b", "a"), 5))
  fit <- markovchainFit(seq, method = "mle", confint = TRUE)
  
  expect_lt(fit$lowerEndpointMatrix["a", "b"], fit$upperEndpointMatrix["a", "b"])
  expect_gt(fit$upperEndpointMatrix["a", "b"], 0)
  expect_lt(fit$lowerEndpointMatrix["a", "a"], 1)
  expect_lte(fit$upperEndpointMatrix["a", "a"], 1)
})

test_that("MLE standard error uses binomial row denominator", {
  seq <- c("a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "b")
  fit <- markovchainFit(seq, method = "mle", confint = TRUE)
  p <- fit$estimate@transitionMatrix["a", "a"]
  n <- sum(createSequenceMatrix(seq, toRowProbs = FALSE)["a", ])
  expect_equal(fit$standardError["a", "a"], sqrt(p * (1 - p) / n), tolerance = 1e-12)
})

test_that("Bootstrap intervals do not shrink deterministically as 1/sqrt(nboot)", {
  set.seed(123)
  seq <- rep(c("a", "b", "a", "a", "b", "b", "a"), 20)
  fit_50 <- markovchainFit(seq, method = "bootstrap", nboot = 50, confint = TRUE)
  set.seed(123)
  fit_500 <- markovchainFit(seq, method = "bootstrap", nboot = 500, confint = TRUE)

  width_50 <- fit_50$confidenceInterval$upperEndpointMatrix["a", "b"] -
    fit_50$confidenceInterval$lowerEndpointMatrix["a", "b"]
  width_500 <- fit_500$confidenceInterval$upperEndpointMatrix["a", "b"] -
    fit_500$confidenceInterval$lowerEndpointMatrix["a", "b"]

  expect_gt(width_500 / width_50, 0.4)
})

test_that("MAP rejects user supplied 1x1 hyperparam instead of treating it as default", {
  seq <- c("a", "b", "a", "c")
  expect_error(markovchainFit(seq, method = "map", hyperparam = matrix(1, 1, 1)))
})

test_that("MAP credible intervals are equal-tailed beta intervals", {
  seq <- c("a", "a", "b", "a", "b", "b", "a")
  hp <- matrix(1, 2, 2, dimnames = list(c("a", "b"), c("a", "b")))
  fit <- markovchainFit(seq, method = "map", hyperparam = hp, confidencelevel = 0.95)
  freq <- createSequenceMatrix(seq, toRowProbs = FALSE)
  p_shape <- freq["a", "a"] + hp["a", "a"]
  q_shape <- sum(freq["a", ]) - freq["a", "a"] + sum(hp["a", ]) - hp["a", "a"]
  expect_equal(fit$confidenceInterval$lowerEndpointMatrix["a", "a"], qbeta(0.025, p_shape, q_shape), tolerance = 1e-8)
  expect_equal(fit$confidenceInterval$upperEndpointMatrix["a", "a"], qbeta(0.975, p_shape, q_shape), tolerance = 1e-8)
})
