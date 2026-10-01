context("fitMTD")

.readExtdata <- function(file, column) {
  path <- system.file("extdata", file, package = "markovchain")
  skip_if(!nzchar(path), paste(file, "not installed"))
  read.csv(path)[[column]]
}

# log-likelihood of an MTD model computed directly from its definition
.mtdLogLik <- function(x, lambda, Q, start) {
  x <- as.character(x)
  sum(vapply(start:length(x), function(t)
    log(sum(lambda * Q[cbind(x[t - seq_along(lambda)], x[t])])), 0))
}

test_that("fitMTD returns a consistent, valid fit", {
  set.seed(123)
  mc <- new("markovchain", states = c("a", "b", "c"),
            transitionMatrix = matrix(c(0.7, 0.2, 0.1,
                                        0.3, 0.4, 0.3,
                                        0.1, 0.3, 0.6), 3, byrow = TRUE))
  x <- rmarkovchain(300, mc, t0 = "a")
  fit <- fitMTD(x, order = 2)
  expect_named(fit$lambda, c("lag1", "lag2"))
  expect_equal(sum(fit$lambda), 1)
  expect_true(all(fit$lambda >= 0))
  expect_s4_class(fit$estimate, "markovchain")
  expect_equal(unname(rowSums(fit$estimate@transitionMatrix)), rep(1, 3))
  expect_true(fit$converged)
  expect_equal(fit$npar, 3 * 2 + 1)
  expect_equal(fit$nobs, 298)
  expect_equal(fit$start, 3L)
  expect_identical(fit$model, "MTD")
  expect_length(fit$Q, 2)
  expect_equal(fit$Q[[1]], t(fit$estimate@transitionMatrix))
  expect_equal(fit$logLikelihood,
               .mtdLogLik(x, fit$lambda, fit$estimate@transitionMatrix, 3))
  expect_equal(fit$BIC, -2 * fit$logLikelihood + log(fit$nobs) * fit$npar)
  expect_equal(fit$AIC, -2 * fit$logLikelihood + 2 * fit$npar)
})

test_that("an MTD model of order one is the first-order Markov chain", {
  set.seed(2)
  x <- sample(c("x", "y", "z"), 200, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  fit <- fitMTD(x, order = 1)
  mle <- markovchainFit(x)$estimate@transitionMatrix
  expect_equal(fit$estimate@transitionMatrix, mle[rownames(fit$estimate@transitionMatrix),
                                                  colnames(fit$estimate@transitionMatrix)],
               tolerance = 1e-8)
  expect_equal(unname(fit$lambda), 1)
})

test_that("fitMTD is consistent with higherOrderLogLik", {
  set.seed(3)
  x <- sample(letters[1:4], 250, replace = TRUE)
  fit <- fitMTD(x, order = 3, start = 5)
  res <- higherOrderLogLik(x, fit, start = 5)
  expect_equal(res$logLik, fit$logLikelihood, tolerance = 1e-10)
  expect_equal(res$npar, fit$npar)
  expect_equal(res$BIC, fit$BIC, tolerance = 1e-10)
})

test_that("fitMTD reaches the maximum found by a general-purpose optimizer", {
  set.seed(4)
  # simulate from an MTD(2) on two states
  Q <- matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE, dimnames = list(c("0", "1"), c("0", "1")))
  lambda <- c(0.6, 0.4)
  x <- c("0", "1")
  for (t in 3:400) {
    p1 <- sum(lambda * Q[x[t - 1:2], "1"])
    x[t] <- if (runif(1) < p1) "1" else "0"
  }
  fit <- fitMTD(x, order = 2, nstart = 3)
  negLogLik <- function(p) {
    l1 <- plogis(p[1]); q <- plogis(p[2:3])
    Qp <- cbind(1 - q, q); dimnames(Qp) <- dimnames(Q)
    -.mtdLogLik(x, c(l1, 1 - l1), Qp, 3)
  }
  best <- min(vapply(1:5, function(i)
    optim(rnorm(3), negLogLik, control = list(reltol = 1e-12, maxit = 5000))$value, 0))
  expect_gte(fit$logLikelihood, -best - 1e-4)
})

test_that("more starting points never give a lower likelihood", {
  set.seed(5)
  x <- sample(c("a", "b", "c"), 150, replace = TRUE)
  one <- fitMTD(x, order = 3)
  set.seed(6)
  many <- fitMTD(x, order = 3, nstart = 5)
  expect_gte(many$logLikelihood, one$logLikelihood - 1e-8)
})

test_that("the row of a state never used as a departure state is uniform", {
  x <- c("a", "b", "a", "b", "a", "a", "b", "c")
  set.seed(7)
  for (nstart in c(1, 3)) {
    fit <- fitMTD(x, order = 2, nstart = nstart)
    expect_equal(unname(fit$estimate@transitionMatrix["c", ]), rep(1 / 3, 3))
  }
})

test_that("fitMTD accepts numeric and factor sequences", {
  x <- c(1, 2, 2, 1, 1, 2, 1, 2, 2, 2, 1, 1, 2, 1)
  expect_equal(fitMTD(x)$logLikelihood, fitMTD(as.character(x))$logLikelihood)
  expect_equal(fitMTD(factor(x))$logLikelihood, fitMTD(as.character(x))$logLikelihood)
})

test_that("fitMTD validates its arguments", {
  x <- c("a", "b", "a", "a", "b", "b", "a")
  expect_error(fitMTD(c("a", NA, "b", "a")), "missing")
  expect_error(fitMTD(c("a", "b")), "at least three")
  expect_error(fitMTD(x, order = 0), "order")
  expect_error(fitMTD(x, order = 1.5), "order")
  expect_error(fitMTD(x, order = 7), "order")
  expect_error(fitMTD(x, order = 2, start = 2), "start")
  expect_error(fitMTD(x, order = 2, start = 8), "start")
  expect_error(fitMTD(x, nstart = 0), "nstart")
  expect_error(fitMTD(x, tol = 0), "tol")
  expect_error(fitMTD(rep("a", 5)), "two different states")
  expect_warning(fitMTD(x, order = 2, maxit = 1), "did not converge")
})

test_that("fitMTD reproduces the MTD rows of Berchtold and Raftery (2002), Tables 2 and 3", {
  # the paper conditions on the first 14 observations (start = 15) and does not
  # count the elements of Q estimated as zero in the number of parameters
  check <- function(x, order, ll, bic) {
    fit <- fitMTD(x, order = order, start = 15)
    zeros <- sum(fit$estimate@transitionMatrix < 1e-8)
    expect_lt(abs(fit$logLikelihood - ll), 0.06)
    expect_lt(abs(-2 * fit$logLikelihood + (fit$npar - zeros) * log(fit$nobs) - bic), 0.06)
    fit
  }
  wind <- .readExtdata("koeberg_wind.csv", "state")
  fit2 <- check(wind, 2, -393.4, 859.3)
  check(wind, 3, -393.2, 865.6)
  seizures <- .readExtdata("epileptic_seizures.csv", "seizure")
  check(seizures, 2, -119.5, 254.7)
  check(seizures, 3, -117.7, 256.4)

  # estimates printed in Section 1.3 of the paper (rows = departure state)
  expect_equal(unname(fit2$lambda), c(0.7569, 0.2431), tolerance = 5e-4)
  Qprinted <- matrix(c(0.8301, 0.0689, 0.0077, 0.0933,
                       0.0369, 0.9012, 0.0619, 0.0000,
                       0.0155, 0.1553, 0.8070, 0.0222,
                       0.0779, 0.0000, 0.0528, 0.8693), 4, 4, byrow = TRUE)
  expect_lt(max(abs(unname(fit2$estimate@transitionMatrix) - Qprinted)), 1e-3)
})
