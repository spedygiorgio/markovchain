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
  skip_if_not_installed("Rsolnp")
  data(rain)
  for (k in 1:3) {
    fit <- fitHigherOrder(rain$rain, order = k)
    expect_length(fit$lambda, k)
    expect_true(all(fit$lambda >= -1e-6))
    expect_equal(sum(fit$lambda), 1, tolerance = 1e-6)
  }
})

test_that("fitHigherOrder minimizes the distance to the stationary distribution", {
  skip_if_not_installed("Rsolnp")
  data(rain)
  data(preproglucacon)
  for (s in list(rain$rain, preproglucacon$preproglucacon)) {
    fit <- fitHigherOrder(s, order = 2)
    grid <- seq(0, 1, by = 0.01)
    gridObjective <- sapply(grid, function(a) .hoObjective(s, c(a, 1 - a)))
    # the fit must be at least as good as every point of a fine grid
    expect_lte(.hoObjective(s, fit$lambda), min(gridObjective) * (1 + 1e-3))
  }
})

test_that("fitHigherOrder no longer stays at the starting point", {
  skip_if_not_installed("Rsolnp")
  data(preproglucacon)
  fit <- fitHigherOrder(preproglucacon$preproglucacon, order = 2)
  # the minimum of the distance is at about (0.69, 0.31), not at (0.5, 0.5)
  expect_gt(abs(fit$lambda[1] - 0.5), 0.1)
  expect_lt(.hoObjective(preproglucacon$preproglucacon, fit$lambda),
            .hoObjective(preproglucacon$preproglucacon, c(0.5, 0.5)))
})

test_that("order 1 gives weight 1", {
  skip_if_not_installed("Rsolnp")
  expect_equal(fitHigherOrder(c("a","b","a","c","b","a","b","c","a"), order = 1)$lambda, 1,
               tolerance = 1e-6)
})
