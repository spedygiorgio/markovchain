context("higherOrderLogLik")

seqA <- c("a", "a", "b", "b", "a", "c", "b", "a", "b", "c", "a", "b",
          "c", "a", "b", "c", "a", "b", "a", "b")

# a hand-built "fit", so that the tests do not depend on Rsolnp
handFit <- function(s, lambda) {
  Q <- lapply(seq_along(lambda), function(o) seq2matHigh(s, o))
  list(lambda = lambda, Q = Q, X = seq2freqProb(s))
}

test_that("the log-likelihood equals the explicit sum over observations", {
  fit <- handFit(seqA, c(0.6, 0.4))
  states <- rownames(fit$Q[[1]])
  manual <- 0
  for (t in 3:length(seqA)) {
    p <- 0.6 * fit$Q[[1]][seqA[t], seqA[t - 1]] + 0.4 * fit$Q[[2]][seqA[t], seqA[t - 2]]
    manual <- manual + log(p)
  }
  res <- higherOrderLogLik(seqA, fit)
  expect_equal(res$logLik, manual)
  expect_equal(res$deviance, -2 * manual)
  expect_equal(res$nobs, length(seqA) - 2)
  expect_equal(res$order, 2)
  expect_equal(res$start, 3L)
})

test_that("order 1 reduces to the first-order maximum likelihood log-likelihood", {
  fit <- handFit(seqA, 1)
  N <- createSequenceMatrix(seqA)         # rows = from, columns = to
  P <- N / rowSums(N)
  expected <- sum(N[N > 0] * log(P[N > 0]))
  expect_equal(higherOrderLogLik(seqA, fit)$logLik, expected)
})

test_that("AIC and BIC follow from the log-likelihood and the parameter count", {
  fit <- handFit(seqA, c(0.5, 0.5))
  res <- higherOrderLogLik(seqA, fit)
  r <- 3
  npar <- 2 * r * (r - 1) + 1
  expect_equal(res$npar, npar)
  expect_equal(res$AIC, -2 * res$logLik + 2 * npar)
  expect_equal(res$BIC, -2 * res$logLik + log(res$nobs) * npar)
})

test_that("start puts models of different orders on the same observations", {
  f1 <- handFit(seqA, 1)
  f2 <- handFit(seqA, c(0.5, 0.5))
  expect_equal(higherOrderLogLik(seqA, f1, start = 3)$nobs,
               higherOrderLogLik(seqA, f2, start = 3)$nobs)
  expect_gt(higherOrderLogLik(seqA, f1)$nobs, higherOrderLogLik(seqA, f2)$nobs)
})

test_that("it works on the output of fitHigherOrder", {
  skip_if_not_installed("Rsolnp")
  fit <- fitHigherOrder(seqA, order = 2)
  res <- higherOrderLogLik(seqA, fit)
  expect_true(is.finite(res$logLik))
  expect_lt(res$logLik, 0)
  # fit = NULL fits the model internally and gives the same answer
  expect_equal(higherOrderLogLik(seqA, order = 2)$logLik, res$logLik)
})

test_that("a probability of zero gives -Inf instead of NaN", {
  fit <- handFit(seqA, 1)
  fit$Q[[1]][] <- 0
  fit$Q[[1]][1, ] <- 1          # only ever moves to the first state
  res <- higherOrderLogLik(seqA, fit)
  expect_identical(res$logLik, -Inf)
  expect_identical(res$deviance, Inf)
})

test_that("invalid input is rejected", {
  fit <- handFit(seqA, c(0.5, 0.5))
  expect_error(higherOrderLogLik(1:5, fit), "character")
  expect_error(higherOrderLogLik(seqA, list(a = 1)), "fitHigherOrder")
  expect_error(higherOrderLogLik(c(seqA, "z"), fit), "unknown")
  expect_error(higherOrderLogLik(seqA, fit, start = 2), "start")
  expect_error(higherOrderLogLik(seqA, fit, start = length(seqA) + 1), "start")
})
