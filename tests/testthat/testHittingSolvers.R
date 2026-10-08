context("hittingProbabilities(): choice of linear solver (#203)")

# Gambler's ruin: states s1..sn, s1 and sn absorbing, p of moving up.
# The absorption probability towards s1 is known in closed form, which gives
# an external reference independent of any of the three solvers.
.ruinChain <- function(n, p) {
  P <- matrix(0, n, n)
  P[1, 1] <- 1
  P[n, n] <- 1
  for (i in 2:(n - 1)) {
    P[i, i + 1] <- p
    P[i, i - 1] <- 1 - p
  }
  sn <- paste0("s", seq_len(n))
  dimnames(P) <- list(sn, sn)
  new("markovchain", states = sn, transitionMatrix = P)
}

# P(absorbed at 0 | start at fortune k), fortune k <-> state s(k+1)
.ruinExact <- function(n, p) {
  N <- n - 1
  ratio <- (1 - p) / p
  k <- 0:N
  if (abs(p - 0.5) < 1e-15) 1 - k / N else (ratio^k - ratio^N) / (1 - ratio^N)
}

test_that("every solver matches the closed-form gambler's ruin probabilities", {
  for (n in c(25L, 120L)) {
    for (p in c(0.3, 0.49, 0.5)) {
      expected <- .ruinExact(n, p)
      mc <- .ruinChain(n, p)
      for (solver in c("direct", "bicgstab", "doubling")) {
        got <- suppressWarnings(
          hittingProbabilities(mc, targets = "s1", solver = solver,
                               tol = 1e-12, maxIter = 50000)[, 1]
        )
        expect_equal(as.numeric(got), expected, tolerance = 1e-8,
                     info = paste("solver", solver, "n", n, "p", p))
      }
    }
  }
})

test_that("the solvers agree with each other on the full matrix", {
  mc <- .ruinChain(40L, 0.49)
  direct <- hittingProbabilities(mc, solver = "direct")
  doubling <- hittingProbabilities(mc, solver = "doubling")
  bicgstab <- suppressWarnings(
    hittingProbabilities(mc, solver = "bicgstab", tol = 1e-12, maxIter = 50000)
  )
  expect_equal(direct, doubling, tolerance = 1e-7)
  expect_equal(direct, bicgstab, tolerance = 1e-7)
  # dimnames and shape must not depend on the solver
  expect_identical(dimnames(direct), dimnames(doubling))
  expect_identical(dimnames(direct), dimnames(bicgstab))
})

test_that("the default solver satisfies the defining equation to machine precision", {
  # h(i, j) = P(i, j) + sum_{k != j} P(i, k) h(k, j) for every i != j.
  # This is an external check: it uses only the transition matrix and the
  # returned matrix, not the internals of any solver.
  mc <- .ruinChain(60L, 0.49)
  H <- hittingProbabilities(mc)
  P <- mc@transitionMatrix
  n <- nrow(P)
  worst <- 0
  for (j in seq_len(n)) {
    Pm <- P
    Pm[, j] <- 0
    residual <- abs(H[, j] - (P[, j] + as.numeric(Pm %*% H[, j])))[-j]
    worst <- max(worst, max(residual))
  }
  expect_lt(worst, 1e-12)
})

test_that("solver choice is validated and the default is 'direct'", {
  mc <- .ruinChain(12L, 0.4)
  expect_error(hittingProbabilities(mc, solver = "gauss"), "arg")
  expect_identical(hittingProbabilities(mc),
                   hittingProbabilities(mc, solver = "direct"))
  # partial matching of the solver name, as match.arg allows
  expect_identical(hittingProbabilities(mc, solver = "bicg"),
                   hittingProbabilities(mc, solver = "bicgstab"))
})

test_that("tol and maxIter are validated", {
  mc <- .ruinChain(12L, 0.4)
  expect_error(hittingProbabilities(mc, tol = 0), "tol")
  expect_error(hittingProbabilities(mc, tol = -1), "tol")
  expect_error(hittingProbabilities(mc, tol = c(1e-5, 1e-6)), "tol")
  expect_error(hittingProbabilities(mc, tol = NA_real_), "tol")
  expect_error(hittingProbabilities(mc, maxIter = 0), "maxIter")
  expect_error(hittingProbabilities(mc, maxIter = 2.5), "maxIter")
  expect_error(hittingProbabilities(mc, maxIter = NA_integer_), "maxIter")
})

test_that("an iterative solver stopped early warns and reports the residual", {
  mc <- .ruinChain(80L, 0.5)
  # One iteration cannot solve an 78-state system: the warning must name the
  # state and report the residual in scientific notation (a plain
  # std::to_string() would print every residual below 1e-6 as "0.000000").
  expect_warning(
    hittingProbabilities(mc, targets = "s1", solver = "bicgstab", maxIter = 1L),
    "did not fully converge"
  )
  warningText <- tryCatch(
    hittingProbabilities(mc, targets = "s1", solver = "doubling",
                         tol = 1e-15, maxIter = 1L),
    warning = function(w) conditionMessage(w)
  )
  expect_match(warningText, "s1")
  expect_match(warningText, "[0-9]e[+-][0-9]+")
})

test_that("solver works with targets, byrow = FALSE and absorbing states", {
  mc <- .ruinChain(30L, 0.45)
  byCol <- new("markovchain", states = mc@states, byrow = FALSE,
               transitionMatrix = t(mc@transitionMatrix))

  for (solver in c("direct", "bicgstab", "doubling")) {
    rowwise <- suppressWarnings(
      hittingProbabilities(mc, targets = c("s1", "s30"), solver = solver,
                           tol = 1e-12, maxIter = 50000)
    )
    colwise <- suppressWarnings(
      hittingProbabilities(byCol, targets = c("s1", "s30"), solver = solver,
                           tol = 1e-12, maxIter = 50000)
    )
    # a column-stochastic chain returns the transpose, as the matrix does
    expect_equal(rowwise, t(colwise), tolerance = 1e-10,
                 info = paste("solver", solver))
    # the two absorbing states keep their boundary values
    expect_equal(unname(rowwise["s1", "s1"]), 1)
    expect_equal(unname(rowwise["s30", "s30"]), 1)
    expect_equal(unname(rowwise["s1", "s30"]), 0)
    expect_equal(unname(rowwise["s30", "s1"]), 0)
    # the two columns add up to one on every transient state: the chain is
    # absorbed at one barrier or the other with probability one
    transient <- paste0("s", 2:29)
    expect_equal(as.numeric(rowwise[transient, "s1"] + rowwise[transient, "s30"]),
                 rep(1, length(transient)), tolerance = 1e-9)
  }
})

test_that("solver choice does not disturb an irreducible chain", {
  # With a single communicating class the graph pass already decides every
  # state, so the solvers must all return the all-ones matrix.
  P <- matrix(c(0.5, 0.3, 0.2,
                0.1, 0.6, 0.3,
                0.2, 0.2, 0.6), nrow = 3, byrow = TRUE,
              dimnames = list(letters[1:3], letters[1:3]))
  mc <- new("markovchain", transitionMatrix = P)
  for (solver in c("direct", "bicgstab", "doubling")) {
    expect_equal(hittingProbabilities(mc, solver = solver),
                 hittingProbabilities(mc, solver = "direct"),
                 info = paste("solver", solver))
  }
})

test_that("bicgstab survives a Lanczos breakdown (five-state vignette chain)", {
  # On this chain (rHat, r) vanishes exactly at the second iteration for
  # target "5": the solver must restart, not return the interrupted iterate.
  M <- markovchain:::zeros(5)
  M[1, 1] <- M[5, 5] <- 1
  M[2, 1] <- M[2, 3] <- 1/2
  M[3, 2] <- M[3, 4] <- 1/2
  M[4, 2] <- M[4, 5] <- 1/2
  mc <- new("markovchain", transitionMatrix = M)
  direct <- hittingProbabilities(mc)
  expect_silent(bicg <- hittingProbabilities(mc, solver = "bicgstab"))
  expect_equal(bicg, direct, tolerance = 1e-12)
  expect_equal(unname(direct[2:4, 5]), c(.2, .4, .6), tolerance = 1e-12)
})

test_that("bicgstab agrees with the direct solver on random absorbing chains", {
  set.seed(285)
  for (rep in 1:60) {
    n <- sample(3:25, 1)
    P <- matrix(runif(n * n) * (runif(n * n) < runif(1, .1, .6)), n)
    ab <- sample(n, sample(1:max(1, n %/% 4), 1))
    P[ab, ] <- 0
    P[cbind(ab, ab)] <- 1
    z <- which(rowSums(P) == 0)
    P[cbind(z, z)] <- 1
    if (rep %% 3 == 0) P <- round(4 * P) + diag(n) * (rowSums(round(4 * P)) == 0)
    P <- P / rowSums(P)
    nm <- paste0("s", seq_len(n))
    dimnames(P) <- list(nm, nm)
    mc <- new("markovchain", transitionMatrix = P)
    expect_silent(bicg <- hittingProbabilities(mc, solver = "bicgstab"))
    expect_equal(bicg, hittingProbabilities(mc), tolerance = 1e-10, info = rep)
  }
})
