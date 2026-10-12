context("Higher order multivariate Markov chains: Ching, Ng and Fung (2008)")

# Section 4.1 of Ching, W.-K., Ng, M. K. and Fung, E. S. (2008). Higher-order
# multivariate Markov chains and their applications. Linear Algebra and its
# Applications, 428(2-3), 492-507: two sequences, four states, order 2.
.chingData <- function() {
  cbind(c(2, 1, 3, 3, 4, 3, 2, 1, 3, 3, 2, 1),
        c(2, 4, 4, 4, 4, 2, 3, 3, 1, 4, 3, 3))
}

# published transition counts F^(jk)_h, rows = next state of sequence j,
# columns = state of sequence k h steps before
.chingCounts <- list(
  "11_1" = rbind(c(0, 3, 0, 0), c(0, 0, 2, 0), c(2, 0, 2, 1), c(0, 0, 1, 0)),
  "11_2" = rbind(c(0, 0, 2, 0), c(0, 0, 1, 1), c(2, 2, 1, 0), c(0, 0, 1, 0)),
  "12_1" = rbind(c(0, 1, 2, 0), c(0, 1, 0, 1), c(1, 0, 1, 3), c(0, 0, 0, 1)),
  "12_2" = rbind(c(0, 1, 0, 1), c(1, 0, 0, 1), c(0, 1, 2, 2), c(0, 0, 0, 1)),
  "21_1" = rbind(c(1, 0, 0, 0), c(0, 0, 0, 1), c(0, 2, 2, 0), c(1, 1, 3, 0)),
  # printed with 0 in row 4, column 3: that table has 9 pairs instead of 10,
  # the pair (state 3 of S1 at time 3, state 4 of S2 at time 5) is missing; the
  # corrected count is used
  "21_2" = rbind(c(0, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 3, 1), c(2, 1, 1, 0)),
  "22_1" = rbind(c(0, 0, 1, 0), c(0, 0, 0, 1), c(0, 1, 2, 1), c(1, 1, 0, 3)),
  "22_2" = rbind(c(0, 0, 1, 0), c(0, 0, 0, 1), c(1, 1, 0, 2), c(0, 1, 1, 2)))

test_that("the transition matrices and frequencies are those of the paper", {
  d <- t(matrix(as.character(.chingData()), ncol = 2))
  P <- markovchain:::.allTransMat(d, order = 2)
  # P is stored as P1(1,1), P2(1,1), P1(1,2), P2(1,2), ..., P2(2,2): index
  # 2 * 2 * (j - 1) + 2 * (k - 1) + h maps sequence k to sequence j
  for (j in 1:2) for (k in 1:2) for (h in 1:2) {
    Fjk <- .chingCounts[[paste0(j, k, "_", h)]]
    expect_equal(unname(P[, , 4 * (j - 1) + 2 * (k - 1) + h]),
                 sweep(Fjk, 2L, colSums(Fjk), "/"),
                 info = paste0("P(", j, k, ")_", h))
  }
  X <- markovchain:::.allFreqProbMat(d)
  expect_equal(unname(X[, 1, 1]), c(1 / 4, 1 / 4, 5 / 12, 1 / 12))
  expect_equal(unname(X[, 1, 2]), c(1 / 12, 1 / 6, 1 / 3, 5 / 12))
})

test_that("numeric states are matched by label, and arguments are checked", {
  skip_if_not_installed("Rsolnp")
  d <- .chingData()
  fitChar <- suppressMessages(fitHighOrderMultivarMC(matrix(as.character(d), ncol = 2), 2))
  # codes 0-3 used to be matched by position, giving wrong frequencies
  fitNum <- suppressMessages(fitHighOrderMultivarMC(d - 1, 2))
  expect_equal(fitNum@states, c("0", "1", "2", "3"))
  expect_equal(unname(fitNum@P), unname(fitChar@P))
  expect_equal(fitNum@Lambda, fitChar@Lambda, tolerance = 1e-6)
  # codes above the number of states used to be an error
  expect_s4_class(suppressMessages(fitHighOrderMultivarMC(d + 10, 2)), "hommc")
  expect_error(fitHighOrderMultivarMC(d, 0), "order")
  expect_error(fitHighOrderMultivarMC(d, 12), "order")
  d[3, 1] <- NA
  expect_error(fitHighOrderMultivarMC(d, 2), "missing")
})
