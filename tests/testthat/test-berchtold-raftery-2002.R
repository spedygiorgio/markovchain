context("Berchtold & Raftery (2002), Tables 2 and 3")

### Published benchmark for higherOrderLogLik() and assessIndependence().
###
### Berchtold, A. and Raftery, A. E. (2002), Statistical Science 17(3),
### 328-356, Tables 2 and 3 report log-likelihood and BIC for independence,
### Markov chains of order 1-3 and MTD models on two series, kindly provided
### by the authors (see inst/extdata/README_berchtold_raftery.md). Their
### convention is to condition every model on the first 14 observations, so
### that all models are evaluated on the same n - 14 observations; the BIC
### uses n - 14, and the number of parameters counts only those not forced to
### zero. The published values are rounded to one decimal.

L <- 14L

.readSeries <- function(file, column) {
  path <- system.file("extdata", file, package = "markovchain")
  skip_if(!nzchar(path), paste(file, "not installed"))
  read.csv(path)[[column]]
}

# log-likelihood and number of free parameters of an order-k Markov chain
# (k = 0: independence) estimated by maximum likelihood on times t = L+1..n
.mcFit <- function(x, k, l = L) {
  n <- length(x)
  t <- (l + 1L):n
  ctx <- if (k == 0L) rep("0", length(t)) else
    vapply(t, function(s) paste(x[(s - k):(s - 1L)], collapse = "-"), "")
  key <- paste(ctx, x[t])
  cnt <- table(key)
  ctxCount <- table(ctx)
  context <- sub(" [^ ]*$", "", names(cnt))
  list(logLik = sum(cnt * log(cnt / as.numeric(ctxCount[context]))),
       npar = sum(tapply(x[t], ctx, function(v) length(unique(v)) - 1L)),
       n = length(t))
}

.bic <- function(fit) -2 * fit$logLik + fit$npar * log(fit$n)

.checkTable <- function(x, published) {
  for (k in 0:3) {
    fit <- .mcFit(x, k)
    expect_equal(fit$npar, published$npar[k + 1L])
    expect_lt(abs(fit$logLik - published$ll[k + 1L]), 0.06)
    expect_lt(abs(.bic(fit) - published$bic[k + 1L]), 0.06)
  }
}

test_that("Markov chain and independence rows of Table 2 (Koeberg wind) are reproduced", {
  x <- .readSeries("koeberg_wind.csv", "state")
  expect_length(x, 744)
  expect_equal(as.integer(table(x)), c(141L, 327L, 144L, 132L))
  .checkTable(x, list(ll = c(-954.8, -413.3, -374.9, -346.2),
                      npar = c(3, 11, 27, 39),
                      bic = c(1929.4, 899.1, 927.9, 949.5)))
})

test_that("Markov chain and independence rows of Table 3 (epileptic seizures) are reproduced", {
  x <- .readSeries("epileptic_seizures.csv", "seizure")
  expect_length(x, 204)
  expect_equal(as.integer(table(x)), c(117L, 87L))
  .checkTable(x, list(ll = c(-129.3, -122.6, -119.3, -117.2),
                      npar = c(1, 2, 4, 8),
                      bic = c(263.9, 255.6, 259.6, 276.3)))
})

test_that("higherOrderLogLik reproduces the published first-order Markov chain log-likelihood", {
  for (case in list(list(file = "koeberg_wind.csv", col = "state", ll = -413.3, bic = 899.1),
                    list(file = "epileptic_seizures.csv", col = "seizure", ll = -122.6, bic = 255.6))) {
    x <- as.character(.readSeries(case$file, case$col))
    # maximum likelihood estimate on the transitions that enter the likelihood
    Q <- seq2matHigh(x[L:length(x)], 1)
    res <- higherOrderLogLik(x, list(lambda = 1, Q = list(Q)), start = L + 1L)
    expect_equal(res$nobs, length(x) - L)
    expect_lt(abs(res$logLik - case$ll), 0.06)
    expect_equal(res$logLik, .mcFit(x, 1)$logLik, tolerance = 1e-10)
  }
})

test_that("higherOrderLogLik reproduces the published MTD(2) log-likelihood and BIC (Koeberg)", {
  x <- as.character(.readSeries("koeberg_wind.csv", "state"))
  # Q and lambda as printed in Section 1.3 of the paper; the printed matrix has
  # the departure state on the rows, whereas fit$Q is stored by column
  Qprinted <- matrix(c(0.8301, 0.0689, 0.0077, 0.0933,
                       0.0369, 0.9012, 0.0619, 0.0000,
                       0.0155, 0.1553, 0.8070, 0.0222,
                       0.0779, 0.0000, 0.0528, 0.8693), 4, 4, byrow = TRUE)
  Q <- t(Qprinted)
  dimnames(Q) <- list(as.character(1:4), as.character(1:4))
  fit <- list(lambda = c(0.7569, 0.2431), Q = list(Q, Q))
  res <- higherOrderLogLik(x, fit, start = L + 1L)

  expect_lt(abs(res$logLik - (-393.4)), 0.1)
  # the paper counts the two structural zeros of Q: 4 * 3 - 2 + (2 - 1) = 11
  nparPublished <- 11
  expect_lt(abs(-2 * res$logLik + nparPublished * log(res$nobs) - 859.3), 0.2)
})

test_that("assessIndependence reproduces the likelihood-ratio statistic implied by the published tables", {
  for (case in list(list(file = "koeberg_wind.csv", col = "state", implied = 2 * (954.8 - 413.3), df = 9),
                    list(file = "epileptic_seizures.csv", col = "seizure", implied = 2 * (129.3 - 122.6), df = 1))) {
    x <- as.character(.readSeries(case$file, case$col))
    # transitions t = L + 1..n, the ones that enter the published likelihoods
    res <- assessIndependence(x[L:length(x)], method = "G", verbose = FALSE)
    expect_equal(unname(res$parameter), case$df)
    # exact identity: G = 2 * (loglik of the first-order chain - loglik of independence)
    expect_equal(unname(res$statistic),
                 2 * (.mcFit(x, 1)$logLik - .mcFit(x, 0)$logLik), tolerance = 1e-10)
    # and agreement with the values implied by the rounded published log-likelihoods
    expect_lt(abs(unname(res$statistic) - case$implied), 0.25)
    expect_lt(res$p.value, 1e-3)
  }
})
