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
