library(testthat)
library(markovchain)

states3 <- c("a", "b", "c")
P3 <- matrix(c(0.5, 0.3, 0.2,
               0.2, 0.6, 0.2,
               0.1, 0.1, 0.8), byrow = TRUE, nrow = 3,
             dimnames = list(states3, states3))
mc3 <- new("markovchain", states = states3, transitionMatrix = P3)

## ---- correctness against finite differences -------------------------------

test_that("sensitivity matches a finite-difference recomputation of pi", {
  S <- sensitivity(mc3, "a")
  eps <- 1e-6

  # Perturb p("a"->"c") up and p("a"->"b") down by eps: a valid zero-sum
  # direction on row "a". The predicted derivative is S["c",] - S["b",].
  P2 <- P3
  P2["a", "c"] <- P2["a", "c"] + eps
  P2["a", "b"] <- P2["a", "b"] - eps
  mc2 <- new("markovchain", states = states3, transitionMatrix = P2)

  numeric_d <- (as.numeric(steadyStates(mc2)) - as.numeric(steadyStates(mc3))) / eps
  predicted_d <- S["c", ] - S["b", ]

  expect_equal(unname(predicted_d), numeric_d, tolerance = 1e-4)
})

test_that("sensitivity holds for a general zero-sum perturbation direction", {
  S <- sensitivity(mc3, "b")
  eps <- 1e-6
  set.seed(123)
  d <- rnorm(3)
  d <- d - mean(d) # enforce sum(d) == 0

  P2 <- P3
  P2["b", ] <- P2["b", ] + eps * d
  mc2 <- new("markovchain", states = states3, transitionMatrix = P2)

  numeric_d <- (as.numeric(steadyStates(mc2)) - as.numeric(steadyStates(mc3))) / eps
  predicted_d <- as.numeric(d %*% S)

  expect_equal(predicted_d, numeric_d, tolerance = 1e-4)
})

test_that("sensitivity accepts a numeric state index equivalent to its name", {
  expect_equal(sensitivity(mc3, "a"), sensitivity(mc3, 1L))
  expect_equal(sensitivity(mc3, "c"), sensitivity(mc3, 3L))
})

test_that("sensitivity scales linearly with pi[state] for a fixed direction", {
  # Scaling pi[k] while keeping the same Z-derived direction should scale S
  # by the same factor: check this directly from the closed form rather
  # than assuming any particular chain's pi is "small enough" to matter,
  # since Z itself can grow arbitrarily large for near-reducible chains
  # (a nearly-decomposable chain can have tiny pi[k] but huge Z entries).
  S_a <- sensitivity(mc3, "a")
  pi <- as.numeric(steadyStates(mc3))
  # S["a", ] / pi["a"] recovers (Z[l=1,] - pi), independent of which state
  # was perturbed; confirm the same relation holds for a different state.
  S_b <- sensitivity(mc3, "b")
  expect_equal(S_a / pi[1], S_b / pi[2], tolerance = 1e-8)
})

test_that("sensitivity supports column-stochastic storage", {
  mc3_col <- new("markovchain", states = states3, byrow = FALSE,
                 transitionMatrix = t(P3))
  expect_equal(sensitivity(mc3_col, "a"), sensitivity(mc3, "a"), tolerance = 1e-10)
})

test_that("sensitivity returns a named n-by-n matrix", {
  S <- sensitivity(mc3, "a")
  expect_equal(dim(S), c(3L, 3L))
  expect_equal(dimnames(S), list(states3, states3))
})

## ---- validation and edge cases ---------------------------------------------

test_that("sensitivity rejects reducible chains and invalid state arguments", {
  reducible <- new("markovchain", states = c("a", "b"),
                   transitionMatrix = matrix(c(1, 0, 0.5, 0.5),
                                              byrow = TRUE, nrow = 2,
                                              dimnames = list(c("a", "b"), c("a", "b"))))
  expect_error(sensitivity(reducible, "a"), "irreducible")

  expect_error(sensitivity(mc3, "z"), "Unknown state")
  expect_error(sensitivity(mc3, 0L), "between 1 and")
  expect_error(sensitivity(mc3, 4L), "between 1 and")
  expect_error(sensitivity(mc3, c("a", "b")), "single state")
})
