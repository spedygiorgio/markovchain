library(testthat)
library(markovchain)

states <- c("a", "b")
P <- matrix(c(0.7, 0.3, 0.1, 0.9),
            byrow = TRUE, nrow = 2,
            dimnames = list(states, states))
mc <- new("markovchain", states = states, transitionMatrix = P)

test_that("slem matches the non-unit eigenvalue of a 2-state chain", {
  # Characteristic polynomial: lambda^2 - 1.6*lambda + 0.6 = 0 -> {1, 0.6}.
  expect_equal(slem(mc), 0.6, tolerance = 1e-12)
})

test_that("slem agrees with a direct eigen() computation, excluding only the unit eigenvalue", {
  lambda <- eigen(P, only.values = TRUE)$values
  unit <- which.min(Mod(lambda - 1))
  expected <- max(Mod(lambda[-unit]))
  expect_equal(slem(mc), expected, tolerance = 1e-12)
})

test_that("spectral gap is one minus slem", {
  expect_equal(spectralGap(mc), 1 - slem(mc), tolerance = 1e-12)
  expect_equal(spectralGap(mc), 0.4, tolerance = 1e-12)
})

test_that("implied timescale of a 2-state chain matches -1/log(slem)", {
  ts <- impliedTimescales(mc)
  expect_named(ts, "tau2")
  expect_equal(unname(ts["tau2"]), -1 / log(slem(mc)), tolerance = 1e-12)
})

test_that("slem, spectral gap and implied timescales support column-stochastic storage", {
  mc_col <- new("markovchain", states = states, byrow = FALSE,
                transitionMatrix = t(P))
  expect_equal(slem(mc_col), slem(mc), tolerance = 1e-12)
  expect_equal(spectralGap(mc_col), spectralGap(mc), tolerance = 1e-12)
  expect_equal(impliedTimescales(mc_col), impliedTimescales(mc), tolerance = 1e-12)
})

test_that("a periodic 2-cycle has slem = 1, spectral gap = 0 and an infinite timescale", {
  periodic <- new("markovchain",
                  states = states,
                  transitionMatrix = matrix(c(0, 1, 1, 0), 2, byrow = TRUE,
                                            dimnames = list(states, states)))
  expect_equal(slem(periodic), 1, tolerance = 1e-12)
  expect_equal(spectralGap(periodic), 0, tolerance = 1e-12)
  ts <- impliedTimescales(periodic)
  expect_named(ts, "tau2")
  expect_equal(unname(ts["tau2"]), Inf)
})

test_that("a 3-cycle with complex eigenvalues is handled through the modulus", {
  states3 <- c("a", "b", "c")
  P3 <- matrix(0, 3, 3, dimnames = list(states3, states3))
  P3["a", "b"] <- 1
  P3["b", "c"] <- 1
  P3["c", "a"] <- 1
  mc3 <- new("markovchain", states = states3, transitionMatrix = P3)

  lambda <- eigen(P3, only.values = TRUE)$values
  expect_true(any(Im(lambda) != 0))
  expect_equal(sort(Mod(lambda), decreasing = TRUE), c(1, 1, 1), tolerance = 1e-12)

  expect_equal(slem(mc3), 1, tolerance = 1e-12)
  expect_equal(spectralGap(mc3), 0, tolerance = 1e-12)
  ts <- impliedTimescales(mc3)
  expect_named(ts, c("tau2", "tau3"))
  expect_true(all(is.infinite(ts)))
})

test_that("slem is zero and spectral gap is one for the one-state chain", {
  one_state <- new("markovchain",
                   states = "a",
                   transitionMatrix = matrix(1, 1, 1,
                                             dimnames = list("a", "a")))
  expect_equal(slem(one_state), 0, tolerance = 1e-12)
  expect_equal(spectralGap(one_state), 1, tolerance = 1e-12)
  ts <- impliedTimescales(one_state)
  expect_length(ts, 0)
  expect_type(ts, "double")
})

test_that("a genuinely absorbing-like near-zero eigenvalue gives a near-zero timescale", {
  # A 3-state chain built so that the non-trivial eigenvalues are small and
  # real, well away from both boundaries.
  states3 <- c("a", "b", "c")
  P3 <- matrix(c(1/3, 1/3, 1/3,
                 1/3, 1/3, 1/3,
                 1/3, 1/3, 1/3), byrow = TRUE, nrow = 3,
               dimnames = list(states3, states3))
  mc3 <- new("markovchain", states = states3, transitionMatrix = P3)
  # This is the iid-uniform chain: non-unit eigenvalues are both exactly 0.
  expect_equal(slem(mc3), 0, tolerance = 1e-12)
  expect_equal(spectralGap(mc3), 1, tolerance = 1e-12)
  ts <- impliedTimescales(mc3)
  expect_true(all(ts == 0))
})

test_that("implied timescales are sorted in decreasing order", {
  states4 <- c("a", "b", "c", "d")
  P4 <- matrix(c(0.85, 0.10, 0.03, 0.02,
                 0.05, 0.80, 0.10, 0.05,
                 0.02, 0.08, 0.70, 0.20,
                 0.01, 0.04, 0.15, 0.80), byrow = TRUE, nrow = 4,
               dimnames = list(states4, states4))
  mc4 <- new("markovchain", states = states4, transitionMatrix = P4)
  ts <- impliedTimescales(mc4)
  expect_named(ts, c("tau2", "tau3", "tau4"))
  expect_true(all(diff(unname(ts)) <= 0))
  expect_equal(unname(ts["tau2"]), -1 / log(slem(mc4)), tolerance = 1e-12)
})

test_that("slem, spectral gap and implied timescales reject reducible chains", {
  reducible <- new("markovchain", states = states,
                   transitionMatrix = matrix(c(1, 0, 0.5, 0.5),
                                              byrow = TRUE, nrow = 2,
                                              dimnames = list(states, states)))
  expect_error(slem(reducible), "irreducible")
  expect_error(spectralGap(reducible), "irreducible")
  expect_error(impliedTimescales(reducible), "irreducible")
})
