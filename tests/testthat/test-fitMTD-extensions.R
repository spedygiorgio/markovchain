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
  expect_error(fitMTD(x, 2, model = "MTDg", negativeWeights = TRUE), "only available")
  expect_error(fitMTD(x, 2, negativeWeights = 1), "negativeWeights")
  expect_error(fitMTD(x, 2, accelerate = NA), "accelerate")
  expect_error(fitMTD(x, 2, se = "yes"), "se must")
  expect_null(fitMTD(x, 2)$se)
})

test_that("MTDg weights are not identified, its transition probabilities are", {
  x <- .mtdSim(1500, 6)
  f <- fitMTD(x, 2, model = "MTDg", se = TRUE)
  # moving mass between the lags leaves every probability unchanged
  A1 <- f$lambda[1] * t(f$Q[[1]])
  A2 <- f$lambda[2] * t(f$Q[[2]])
  shift <- 0.02 * matrix(1 / 3, 3, 3)
  l <- c(sum(A1[1, ]) + 0.02, sum(A2[1, ]) - 0.02)
  g <- list(lambda = l, Q = list(t((A1 + shift) / l[1]), t((A2 - shift) / l[2])))
  expect_equal(higherOrderLogLik(x, g, start = 3)$logLik, f$logLikelihood, tolerance = 1e-10)
  expect_true(all(is.na(f$se$lambda)))
  expect_true(all(f$se$probabilities > 0))
  expect_equal(colnames(f$se$probabilities), c("a", "b", "c"))
  expect_true("a-b" %in% rownames(f$se$probabilities))
})

test_that("standard errors of the transition probabilities agree with the Monte Carlo spread", {
  skip_on_cran()
  st <- letters[1:3]
  Q1 <- matrix(c(.7, .2, .1, .1, .8, .1, .2, .2, .6), 3, byrow = TRUE, dimnames = list(st, st))
  Q2 <- matrix(c(.3, .4, .3, .5, .2, .3, .1, .3, .6), 3, byrow = TRUE, dimnames = list(st, st))
  trg <- list(lambda = c(.6, .4), Q = list(t(Q1), t(Q2)))
  set.seed(10)
  est <- t(vapply(seq_len(100), function(i) {
    x <- higherOrderSimulate(1500, trg, t0 = c("a", "b"))
    f <- fitMTD(x, 2, model = "MTDg", se = TRUE)
    p <- higherOrderPredict(f, rbind(c("a", "b"), c("c", "c")))
    c(p[1, 1], f$se$probabilities["a-b", "a"], p[2, 3], f$se$probabilities["c-c", "c"])
  }, numeric(4)))
  expect_equal(mean(est[, 2]) / stats::sd(est[, 1]), 1, tolerance = 0.25)
  expect_equal(mean(est[, 4]) / stats::sd(est[, 3]), 1, tolerance = 0.25)
})

test_that("negative weights: Raftery's condition, nesting and the likelihood", {
  x <- read.csv(system.file("extdata", "epileptic_seizures.csv",
                            package = "markovchain"))$seizure
  f <- fitMTD(x, 5, start = 15, negativeWeights = TRUE)
  g <- fitMTD(x, 5, start = 15)
  expect_true(f$negativeWeights)
  expect_false(g$negativeWeights)
  expect_gte(f$logLikelihood, g$logLikelihood - 1e-8)
  expect_lt(min(f$lambda), 0)
  expect_equal(sum(f$lambda), 1)
  expect_equal(f$npar, g$npar)
  # every transition probability, for all 2^5 histories, lies in [0, 1]
  hist <- as.matrix(expand.grid(rep(list(c("0", "1")), 5), stringsAsFactors = FALSE))
  p <- higherOrderPredict(f, hist)
  expect_true(all(p >= -1e-12 & p <= 1 + 1e-12))
  # higherOrderPredict() normalizes, so check that the raw mixture sums to one
  raw <- apply(hist, 1L, function(h) sum(vapply(1:5, function(i)
    f$lambda[i] * f$Q[[i]][, h[5 - i + 1]], numeric(2))))
  expect_equal(unname(raw), rep(1, nrow(hist)), tolerance = 1e-10)
  expect_equal(higherOrderLogLik(x, f, start = 15)$logLik, f$logLikelihood, tolerance = 1e-10)
})

test_that("standard errors with negative weights agree with the Monte Carlo spread", {
  skip_on_cran()
  st <- letters[1:3]
  Qn <- matrix(1 / 6, 3, 3, dimnames = list(st, st))
  diag(Qn) <- 2 / 3
  trn <- list(lambda = c(1.2, -0.2), Q = list(t(Qn), t(Qn)))
  set.seed(11)
  est <- t(vapply(seq_len(60), function(i) {
    x <- higherOrderSimulate(1500, trn, t0 = c("a", "b"))
    f <- suppressWarnings(fitMTD(x, 2, negativeWeights = TRUE, se = TRUE))
    c(f$lambda[2], f$se$lambda[2])
  }, numeric(2)))
  expect_equal(mean(est[, 1]), -0.2, tolerance = 0.05, scale = 1)
  expect_equal(mean(est[, 2]) / stats::sd(est[, 1]), 1, tolerance = 0.3)
})
