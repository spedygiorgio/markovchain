library(testthat)
library(markovchain)

states2 <- c("a", "b")
states3 <- c("a", "b", "c")

P2 <- matrix(c(0.7, 0.3, 0.1, 0.9),
             byrow = TRUE, nrow = 2,
             dimnames = list(states2, states2))
mc2 <- new("markovchain", states = states2, transitionMatrix = P2)

## ---- lazyChain -------------------------------------------------------------

test_that("lazyChain matches alpha*I + (1-alpha)*P directly", {
  lazyMc <- lazyChain(mc2, alpha = 0.4)
  expected <- 0.4 * diag(2) + 0.6 * P2
  dimnames(expected) <- dimnames(P2)
  expect_equal(as.matrix(lazyMc@transitionMatrix), expected, tolerance = 1e-12)
})

test_that("lazyChain preserves the stationary distribution", {
  lazyMc <- lazyChain(mc2, alpha = 0.7)
  expect_equal(as.numeric(steadyStates(mc2)), as.numeric(steadyStates(lazyMc)),
               tolerance = 1e-10)
})

test_that("alpha = 0 returns the original chain and alpha = 1 never moves", {
  lazyZero <- lazyChain(mc2, alpha = 0)
  expect_equal(as.matrix(lazyZero@transitionMatrix), P2, tolerance = 1e-12)

  lazyOne <- lazyChain(mc2, alpha = 1)
  expect_equal(unname(as.matrix(lazyOne@transitionMatrix)), diag(2),
               tolerance = 1e-12)
})

test_that("lazyChain removes periodicity from a periodic chain", {
  periodic <- new("markovchain", states = states2,
                  transitionMatrix = matrix(c(0, 1, 1, 0), 2, byrow = TRUE,
                                            dimnames = list(states2, states2)))
  expect_equal(period(periodic), 2L)
  lazyPeriodic <- lazyChain(periodic, alpha = 0.5)
  expect_equal(period(lazyPeriodic), 1L)
})

test_that("lazyChain preserves the byrow storage convention", {
  mc2_col <- new("markovchain", states = states2, byrow = FALSE,
                 transitionMatrix = t(P2))
  lazyCol <- lazyChain(mc2_col, alpha = 0.3)
  expect_false(lazyCol@byrow)
  # Converting back to row-stochastic should match the row-stochastic case.
  lazyRow <- lazyChain(mc2, alpha = 0.3)
  expect_equal(t(as.matrix(lazyCol@transitionMatrix)),
               as.matrix(lazyRow@transitionMatrix), tolerance = 1e-12)
})

test_that("non-trivial eigenvalues of the lazy chain follow alpha + (1-alpha)*lambda", {
  lambda <- eigen(P2, only.values = TRUE)$values
  unit <- which.min(Mod(lambda - 1))
  nonTrivial <- lambda[-unit]
  expectedLazy <- 0.4 + 0.6 * nonTrivial

  lazyMc <- lazyChain(mc2, alpha = 0.4)
  lazyLambda <- eigen(as.matrix(lazyMc@transitionMatrix), only.values = TRUE)$values
  lazyUnit <- which.min(Mod(lazyLambda - 1))
  expect_equal(sort(Mod(lazyLambda[-lazyUnit])), sort(Mod(expectedLazy)),
               tolerance = 1e-10)
})

test_that("lazyChain validates alpha", {
  expect_error(lazyChain(mc2, alpha = -0.1), "alpha")
  expect_error(lazyChain(mc2, alpha = 1.1), "alpha")
  expect_error(lazyChain(mc2, alpha = c(0.2, 0.3)), "alpha")
  expect_error(lazyChain(mc2, alpha = NA_real_), "alpha")
})

## ---- subchain ---------------------------------------------------------------

P3 <- matrix(c(0.5, 0.3, 0.2,
               0.2, 0.6, 0.2,
               0.1, 0.1, 0.8), byrow = TRUE, nrow = 3,
             dimnames = list(states3, states3))
mc3 <- new("markovchain", states = states3, transitionMatrix = P3)

test_that("subchain submatrix is the raw principal submatrix, generally sub-stochastic", {
  sub <- subchain(mc3, c("a", "b"), method = "submatrix")
  expect_true(is.matrix(sub))
  expect_false(is(sub, "markovchain"))
  expect_equal(sub, P3[c("a", "b"), c("a", "b")], tolerance = 1e-12)
  expect_true(all(rowSums(sub) < 1))
})

test_that("subchain renormalize returns a valid markovchain with rows summing to one", {
  watched <- subchain(mc3, c("a", "b"), method = "renormalize")
  expect_true(is(watched, "markovchain"))
  expect_equal(unname(rowSums(as.matrix(watched@transitionMatrix))), c(1, 1),
               tolerance = 1e-12)

  expected <- P3[c("a", "b"), c("a", "b")] / rowSums(P3[c("a", "b"), c("a", "b")])
  expect_equal(unname(as.matrix(watched@transitionMatrix)), unname(expected),
               tolerance = 1e-12)
})

test_that("subchain renormalize errors clearly when a state cannot stay inside the subset", {
  absStates <- c("a", "b", "c")
  absMat <- matrix(c(1, 0, 0,
                     0.5, 0, 0.5,
                     0, 0, 1), byrow = TRUE, nrow = 3,
                   dimnames = list(absStates, absStates))
  absMc <- new("markovchain", states = absStates, transitionMatrix = absMat)
  expect_error(subchain(absMc, "b", method = "renormalize"), "zero probability")
})

test_that("subchain works with a single retained state", {
  single <- subchain(mc3, "a", method = "submatrix")
  expect_equal(dim(single), c(1L, 1L))
  expect_equal(single[1, 1], P3["a", "a"])
})

test_that("subchain supports column-stochastic storage", {
  mc3_col <- new("markovchain", states = states3, byrow = FALSE,
                 transitionMatrix = t(P3))
  subCol <- subchain(mc3_col, c("a", "b"), method = "submatrix")
  subRow <- subchain(mc3, c("a", "b"), method = "submatrix")
  expect_equal(subCol, subRow, tolerance = 1e-12)
})

test_that("subchain validates its states argument", {
  expect_error(subchain(mc3, character(0)), "non-empty")
  expect_error(subchain(mc3, c("a", "a")), "duplicate")
  expect_error(subchain(mc3, c("a", "z")), "Unknown state")
  expect_error(subchain(mc3, NA_character_), "missing")
})

test_that("subchain and lazyChain compose: lazifying a subchain still preserves its stationary distribution", {
  watched <- subchain(mc3, c("a", "b"), method = "renormalize")
  lazyWatched <- lazyChain(watched, alpha = 0.5)
  expect_equal(as.numeric(steadyStates(watched)), as.numeric(steadyStates(lazyWatched)),
               tolerance = 1e-10)
})
