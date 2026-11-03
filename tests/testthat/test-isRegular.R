# is.regular() is decided as "irreducible and aperiodic" (src/probabilistic.cpp); here it is
# checked against the definition (some power of P has only positive entries), computed exactly on
# the 0/1 pattern up to Wielandt's bound m^2 - 2m + 2.

.regularByDefinition <- function(P) {
  Z <- (P > 0) * 1
  m <- nrow(Z)
  A <- Z
  for (k in seq_len(m * m - 2 * m + 2)) {
    if (all(A > 0)) return(TRUE)
    A <- ((A %*% Z) > 0) * 1
  }
  all(A > 0)
}

test_that("is.regular() on textbook cases", {
  mk <- function(P, byrow = TRUE) new("markovchain", transitionMatrix = P, byrow = byrow,
                                       states = paste0("s", seq_len(nrow(P))))
  expect_true(is.regular(mk(matrix(1, 1, 1))))
  expect_true(is.regular(mk(matrix(c(0.5, 0.5, 1, 0), 2, byrow = TRUE))))
  expect_false(is.regular(mk(matrix(c(0, 1, 1, 0), 2, byrow = TRUE))))           # period 2
  expect_false(is.regular(mk(matrix(c(1, 0, 0.5, 0.5), 2, byrow = TRUE))))       # reducible
  cyc <- matrix(c(0, 1, 0, 0, 0, 1, 1, 0, 0), 3, byrow = TRUE)
  expect_false(is.regular(mk(cyc)))
  cyc[1, 1] <- 0.5; cyc[1, 2] <- 0.5                                             # a self-loop
  expect_true(is.regular(mk(cyc)))
  expect_true(is.regular(mk(t(cyc), byrow = FALSE)))
  # Wielandt's extremal matrix: regular, but only its (m^2 - 2m + 2)-th power is positive
  m <- 6
  W <- matrix(0, m, m)
  W[cbind(1:(m - 1), 2:m)] <- 1
  W[m, 1:2] <- 0.5
  expect_true(is.regular(mk(W)))
})

test_that("is.regular() agrees with the definition on random chains", {
  set.seed(42)
  for (s in 1:400) {
    n <- sample(1:8, 1)
    Z <- matrix(runif(n * n) < runif(1, 0.1, 0.6), n)
    if (s %% 3 == 0) {                       # block-cyclic patterns, often periodic
      d <- sample(2:3, 1); g <- sample(rep_len(1:d, n))
      Z <- outer(g, g, function(a, b) b == (a %% d) + 1) & Z
    }
    for (i in which(rowSums(Z) == 0)) Z[i, sample(n, 1)] <- TRUE
    P <- Z * matrix(rexp(n * n), n)
    P <- P / rowSums(P)
    mc <- new("markovchain", transitionMatrix = P, states = paste0("s", seq_len(n)))
    expect_identical(is.regular(mc), .regularByDefinition(P), info = paste("seed step", s))
  }
})

test_that("period() is right on a large periodic chain", {
  # bipartite chain on 400 states: period 2, not regular
  n <- 400
  set.seed(7)
  P <- matrix(0, n, n)
  ev <- seq(1, n, 2); od <- seq(2, n, 2)
  for (i in ev) P[i, sample(od, 3)] <- 1
  for (i in od) P[i, sample(ev, 3)] <- 1
  P[cbind(ev, od)] <- 1
  P[cbind(od[-length(od)], ev[-1])] <- 1
  P <- P / rowSums(P)
  mc <- new("markovchain", transitionMatrix = P, states = paste0("s", seq_len(n)))
  expect_true(is.irreducible(mc))
  expect_identical(period(mc), 2L)
  expect_false(is.regular(mc))
})
