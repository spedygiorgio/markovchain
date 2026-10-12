context("fitMTD: MTDg model, SQUAREM acceleration, standard errors")

.mtdSim <- function(n, seed) {
  set.seed(seed)
  Q <- matrix(c(.7, .2, .1, .1, .8, .1, .2, .2, .6), 3, byrow = TRUE,
              dimnames = list(letters[1:3], letters[1:3]))
  tr <- list(lambda = c(.6, .4), Q = list(t(Q), t(Q)))
  higherOrderSimulate(n, tr, t0 = c("a", "b"))
}

test_that("MTDg dominates MTD and the fitHigherOrder mle model, with consistent npar", {
  x <- .mtdSim(1500, 1)
  g <- fitMTD(x, 3, model = "MTDg")
  m <- fitMTD(x, 3)
  h <- fitHigherOrder(x, 3, method = "mle")
  expect_equal(g$model, "MTDg")
  expect_equal(g$npar, (3 - 1) * (1 + 3 * (3 - 1)))
  expect_gte(g$logLikelihood, m$logLikelihood - 1e-8)
  hl <- higherOrderLogLik(x, h, start = g$start)
  expect_gte(g$logLikelihood, hl$logLik - 1e-6)
  gl <- higherOrderLogLik(x, g, start = g$start)
  expect_equal(gl$logLik, g$logLikelihood, tolerance = 1e-8)
  expect_equal(gl$npar, g$npar)
  expect_length(g$estimate, 3)
  expect_equal(sum(g$lambda), 1)
  for (q in g$Q) expect_equal(unname(colSums(q)), rep(1, 3))
  p <- higherOrderPredict(g, c("a", "b", "c"))
  expect_equal(sum(p), 1)
})

test_that("MTDg likelihood is non-decreasing with the order and matches a direct maximization", {
  x <- .mtdSim(400, 2)
  l <- vapply(1:3, function(k) fitMTD(x, k, start = 4, model = "MTDg")$logLikelihood, numeric(1))
  expect_true(all(diff(l) >= -1e-8))
  # order 1 is the Markov chain of order 1: the empirical likelihood
  mc <- selectOrder(x, maxOrder = 1, start = 4)$table$logLik[2]
  expect_equal(l[1], mc, tolerance = 1e-8)
})

test_that("accelerated EM reaches the same maximum with far fewer steps", {
  x <- .mtdSim(2000, 3)
  for (model in c("MTD", "MTDg")) {
    a <- fitMTD(x, 3, model = model)
    b <- fitMTD(x, 3, model = model, accelerate = TRUE)
    expect_equal(b$logLikelihood, a$logLikelihood, tolerance = 1e-6)
    expect_true(b$converged)
  }
  set.seed(4)
  y <- sample(letters[1:3], 20000, TRUE)
  b <- expect_warning(fitMTD(y, 3), "did not converge")
  a <- fitMTD(y, 3, accelerate = TRUE)
  expect_true(a$converged)
  expect_lt(a$iterations, 5000)
  expect_gte(a$logLikelihood, b$logLikelihood - 1e-4)
})

test_that("standard errors agree with the Monte Carlo spread", {
  skip_on_cran()
  est <- t(vapply(seq_len(100), function(i) {
    f <- fitMTD(.mtdSim(1500, 100 + i), 2, se = TRUE)
    c(f$lambda[1], f$se$lambda[1], f$estimate@transitionMatrix[1, 1], f$se$Q[1, 1])
  }, numeric(4)))
  expect_equal(mean(est[, 2]) / stats::sd(est[, 1]), 1, tolerance = 0.25)
  expect_equal(mean(est[, 4]) / stats::sd(est[, 3]), 1, tolerance = 0.25)
})

test_that("standard errors flag the boundary and arguments are validated", {
  x <- .mtdSim(1500, 5)
  f <- suppressWarnings(fitMTD(x, 3, se = TRUE))
  expect_equal(dim(f$se$Q), c(3, 3))
  expect_error(fitMTD(x, 2, model = "MTDg", se = TRUE), "only available")
  expect_error(fitMTD(x, 2, accelerate = NA), "accelerate")
  expect_error(fitMTD(x, 2, se = "yes"), "se must")
  expect_null(fitMTD(x, 2)$se)
})
