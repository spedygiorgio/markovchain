context("fitHigherOrder / seq2matHigh")

### Regression test for the fix in seq2matHigh() (src/fitHigherOrder.cpp):
### a state that never occurs as the "from" state at a given lag produced
### an entire column of NaN (0/0), which then silently propagated through
### fitHigherOrder()'s Q %*% X products, breaking the fit with no clear
### indication why (a single NaN column contaminates the entire matrix
### product). The fix falls back to a uniform distribution for that
### column instead.

test_that("seq2matHigh reproduces the standard lag-order counts when every state has data", {
  s <- c("a", "b", "c", "a", "b", "c", "a", "b", "c", "a")
  Q <- seq2matHigh(s, 2)

  expect_equal(unname(colSums(Q)), c(1, 1, 1))
  expect_false(any(is.nan(Q)))
})

test_that("seq2matHigh falls back to a uniform column instead of NaN when a state never occurs as 'from' at that lag", {
  # "c" occurs only in the last position, so it is never a "from" state at
  # lag 1 -- before the fix, colsums["c"] == 0 produced 0/0 == NaN for the
  # entire "c" column.
  s <- c("a", "b", "a", "b", "a", "b", "a", "b", "c")
  Q <- seq2matHigh(s, 1)

  expect_false(any(is.nan(Q)))
  expect_equal(unname(Q[, "c"]), rep(1 / 3, 3))
  expect_equal(unname(colSums(Q)), c(1, 1, 1))
})

test_that("the NaN no longer propagates into fitHigherOrder()'s Q %*% X products", {
  # Reproduces the exact failure mode: previously, a single NaN column in
  # Q turned the ENTIRE Q %*% X product into NaN (matrix multiplication
  # mixes every row), which would make fitHigherOrder()'s quadratic
  # program fail silently.
  s <- c("a", "b", "a", "b", "a", "b", "a", "b", "c")
  X <- seq2freqProb(s)

  for (lag in 1:2) {
    Q <- seq2matHigh(s, lag)
    QX <- Q %*% X
    expect_false(any(is.nan(QX)))
  }
})

### Regression tests for the objective function of fitHigherOrder().
### The objective used to be sum_i(lambda_i * Q_i X - X), i.e. it subtracted
### `order` times the stationary distribution X instead of once. That made it
### (almost) constant in lambda, so the optimizer stayed at its starting point:
### for order 2 and 3 the returned weights were exactly 1/2 and 1/3 each. Even
### with the objective fixed, its value is of the order of 1e-7 on real
### sequences, below solnp's default tolerance, so it is also scaled.

.hoObjective <- function(s, lambda) {
  X <- seq2freqProb(s)
  fitted <- 0
  for (i in seq_along(lambda)) fitted <- fitted + lambda[i] * (seq2matHigh(s, i) %*% X)
  sum((fitted - X)^2)
}

test_that("fitHigherOrder returns valid weights", {
  data(rain)
  for (method in c("mle", "lsq")) for (k in 1:3) {
    if (method == "lsq") skip_if_not_installed("Rsolnp")
    fit <- fitHigherOrder(rain$rain, order = k, method = method)
    expect_length(fit$lambda, k)
    expect_true(all(fit$lambda >= -1e-6))
    expect_equal(sum(fit$lambda), 1, tolerance = 1e-6)
  }
  # the default is the maximum likelihood fit
  expect_equal(fitHigherOrder(rain$rain, 2), fitHigherOrder(rain$rain, 2, method = "mle"))
})

test_that("fitHigherOrder minimizes the distance to the stationary distribution", {
  skip_if_not_installed("Rsolnp")
  data(rain)
  data(preproglucacon)
  for (s in list(rain$rain, preproglucacon$preproglucacon)) {
    fit <- fitHigherOrder(s, order = 2, method = "lsq")
    grid <- seq(0, 1, by = 0.01)
    gridObjective <- sapply(grid, function(a) .hoObjective(s, c(a, 1 - a)))
    # the fit must be at least as good as every point of a fine grid
    expect_lte(.hoObjective(s, fit$lambda), min(gridObjective) * (1 + 1e-3))
  }
})

test_that("fitHigherOrder no longer stays at the starting point", {
  skip_if_not_installed("Rsolnp")
  data(preproglucacon)
  fit <- fitHigherOrder(preproglucacon$preproglucacon, order = 2, method = "lsq")
  # the minimum of the distance is at about (0.69, 0.31), not at (0.5, 0.5)
  expect_gt(abs(fit$lambda[1] - 0.5), 0.1)
  expect_lt(.hoObjective(preproglucacon$preproglucacon, fit$lambda),
            .hoObjective(preproglucacon$preproglucacon, c(0.5, 0.5)))
})

test_that("order 1 gives weight 1", {
  skip_if_not_installed("Rsolnp")
  for (method in c("lsq", "mle")) {
    expect_equal(fitHigherOrder(c("a","b","a","c","b","a","b","c","a"), order = 1,
                                method = method)$lambda, 1, tolerance = 1e-6)
  }
})

### method = "mle": weights by maximum likelihood (EM), same lag matrices.

.hoLogLik <- function(s, lambda, Q) higherOrderLogLik(s, list(lambda = lambda, Q = Q))$logLik

test_that("method = 'mle' returns weights on the simplex and the usual structure", {
  data(rain)
  for (k in 1:3) {
    fit <- fitHigherOrder(rain$rain, order = k, method = "mle")
    expect_named(fit, c("lambda", "Q", "X"))
    expect_length(fit$lambda, k)
    expect_true(all(fit$lambda >= 0))
    expect_equal(sum(fit$lambda), 1, tolerance = 1e-10)
    expect_length(fit$Q, k)
  }
  # order 1 has a single weight
  expect_equal(fitHigherOrder(rain$rain, order = 1, method = "mle")$lambda, 1)
})

test_that("method = 'mle' uses the same lag matrices as method = 'lsq'", {
  skip_if_not_installed("Rsolnp")
  data(preproglucacon)
  s <- preproglucacon$preproglucacon
  lsq <- fitHigherOrder(s, 2, method = "lsq")
  mle <- fitHigherOrder(s, 2, method = "mle")
  expect_equal(mle$Q, lsq$Q)
  expect_equal(mle$X, lsq$X)
})

test_that("the default method is the maximum likelihood one", {
  data(rain)
  expect_identical(fitHigherOrder(rain$rain, 2), fitHigherOrder(rain$rain, 2, method = "mle"))
  expect_error(fitHigherOrder(rain$rain, 2, method = "foo"), "should be one of")
})

test_that("method = 'mle' attains the global maximum of the log-likelihood", {
  data(rain)
  data(preproglucacon)
  for (s in list(rain$rain, preproglucacon$preproglucacon)) {
    # order 2: fine grid on the segment
    fit2 <- fitHigherOrder(s, 2, method = "mle")
    grid <- seq(0, 1, by = 0.001)
    gridLL <- vapply(grid, function(a) .hoLogLik(s, c(a, 1 - a), fit2$Q), 0)
    expect_gte(.hoLogLik(s, fit2$lambda, fit2$Q), max(gridLL) - 1e-6)
    expect_lt(abs(fit2$lambda[1] - grid[which.max(gridLL)]), 0.005)

    # order 3: grid on the simplex
    fit3 <- fitHigherOrder(s, 3, method = "mle")
    best <- -Inf
    for (a in seq(0, 1, by = 0.05)) for (b in seq(0, 1 - a, by = 0.05))
      best <- max(best, .hoLogLik(s, c(a, b, max(0, 1 - a - b)), fit3$Q))
    expect_gte(.hoLogLik(s, fit3$lambda, fit3$Q), best - 1e-6)
  }
})

test_that("maximum likelihood weights never give a lower log-likelihood than least squares", {
  skip_if_not_installed("Rsolnp")
  data(rain)
  data(preproglucacon)
  for (s in list(rain$rain, preproglucacon$preproglucacon)) for (k in 2:3) {
    lsq <- fitHigherOrder(s, k, method = "lsq")
    mle <- fitHigherOrder(s, k, method = "mle")
    expect_gte(higherOrderLogLik(s, mle)$logLik, higherOrderLogLik(s, lsq)$logLik - 1e-8)
  }
})

test_that("the weights satisfy the optimality conditions and beat a long EM run", {
  # Optimality (KKT) conditions of the concave problem on the simplex:
  # g_o <= 1 for every lag, with g_o = 1 where lambda_o > 0, where g_o is the
  # derivative of the log-likelihood in lambda_o divided by the number of
  # observations. Checked without grouping the observations, together with
  # the log-likelihood reached by 20000 plain EM iterations (the method used
  # before), which the weights must not lose to.
  check <- function(seqs, order) {
    fit <- fitHigherOrder(seqs, order)
    states <- rownames(fit$Q[[1L]])
    q <- do.call(rbind, lapply(seqs[lengths(seqs) > order], function(s) {
      idx <- match(s, states)
      times <- (order + 1L):length(idx)
      sapply(seq_len(order), function(o) fit$Q[[o]][cbind(idx[times], idx[times - o])])
    }))
    g <- colMeans(q / drop(q %*% fit$lambda))
    info <- paste("order", order)
    expect_true(all(g <= 1 + 1e-8), info = info)
    expect_true(all(abs(g[fit$lambda > 0] - 1) <= 1e-8), info = info)
    expect_equal(sum(fit$lambda), 1)
    lambda <- rep(1 / order, order)
    for (i in seq_len(20000L)) {
      w <- sweep(q, 2L, lambda, "*")
      lambda <- colMeans(w / rowSums(w))
    }
    expect_gte(sum(log(q %*% fit$lambda)), sum(log(q %*% lambda)) - 1e-9)
    fit$lambda
  }
  data(rain, package = "markovchain")
  for (order in 2:5) check(list(rain$rain), order)
  check(list(rain$rain[1:500], rain$rain[501:1096]), 3)
  # i.i.d. data: the lags carry the same information, EM needed more than
  # 10000 iterations at order 4 and stopped short of the maximum
  set.seed(2)
  x <- sample(letters[1:4], 3000, TRUE)
  expect_silent(check(list(x), 4))
  # a first-order chain fitted at order 3: the higher lags are worth nothing,
  # their weights are (essentially) 0 and their gradients at most 1
  set.seed(3)
  s <- c("a", "b", "c")
  mc <- new("markovchain", states = s, transitionMatrix = matrix(
    c(0.8, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8), 3, byrow = TRUE, dimnames = list(s, s)))
  lambda <- check(list(rmarkovchain(3000, mc, t0 = "a")), 3)
  expect_gt(lambda[1], 0.9)
})

test_that("the Newton solver handles boundary, collinear and degenerate problems", {
  newton <- markovchain:::.mixtureWeightsNewton
  kkt <- function(q, counts, l) {
    g <- colSums(q * (counts / drop(q %*% l))) / sum(counts)
    max(g - 1, abs(g[l > 0] - 1))
  }
  set.seed(10)
  for (rep in 1:200) {
    k <- sample(1:8, 1)
    H <- sample(c(3, 5, 20, 200), 1)
    q <- matrix(runif(H * k)^sample(c(1, 3, 8), 1), H, k)
    if (k > 2 && rep %% 3 == 0) q[, 2] <- q[, 1] * (1 + 1e-9 * runif(H))  # nearly collinear
    if (k > 1 && rep %% 4 == 0) q[, k] <- q[, k] * 1e-3                    # a useless lag
    q <- pmax(q, 1e-12)
    counts <- sample(1:50, H, TRUE)
    l <- newton(q, counts)
    expect_true(attr(l, "converged"), info = paste("case", rep))
    expect_lte(kkt(q, counts, as.numeric(l)), 1e-8)
    expect_true(all(l >= 0) && abs(sum(l) - 1) < 1e-12)
  }
  # one component: weight 1
  expect_equal(as.numeric(newton(matrix(c(0.2, 0.5), 2, 1), c(3, 4))), 1)
})
