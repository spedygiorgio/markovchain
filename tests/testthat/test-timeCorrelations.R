context("timeCorrelations() and timeRelaxations()")

# Chain with complex eigenvalues (-0.35 +/- 0.606i). Reference values from
# numpy, sum(f * pi * (P^t g)) and mu P^t f with numpy.linalg.matrix_power.
# PyDTMC 9.0.0 gets this case wrong as soon as a time point exceeds the
# number of states: it then uses an eigendecomposition whose left and right
# eigenvectors are not biorthonormal and returns 4.0 at every lag.
P <- matrix(c(0.1, 0.8, 0.1,
              0.1, 0.1, 0.8,
              0.8, 0.1, 0.1), 3, byrow = TRUE,
            dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
mc <- new("markovchain", transitionMatrix = P)
x <- c("a", "b", "b", "c", "a", "a")
y <- c("c", "c", "b")
tp <- c(0, 1, 2, 5, 10)

test_that("autocorrelation and cross-correlation match the reference", {
  expect_equal(unname(timeCorrelations(mc, x, timePoints = tp)),
               c(4.666666666666667, 3.766666666666667, 3.836666666666668,
                 3.943976666666668, 3.990584158366669), tolerance = 1e-12)
  expect_equal(unname(timeCorrelations(mc, x, y, timePoints = tp)),
               c(1.3333333333333335, 2.233333333333334, 2.1633333333333336,
                 2.056023333333334, 2.0094158416333343), tolerance = 1e-12)
  expect_named(timeCorrelations(mc, x, timePoints = tp), as.character(tp))
})

test_that("relaxations match the reference", {
  expect_equal(unname(timeRelaxations(mc, x, initial = c(0.5, 0.25, 0.25),
                                      timePoints = tp)),
               c(2.25, 2.0, 1.8775000000000002, 1.9579825000000004,
                 2.000000000000001), tolerance = 1e-12)
  # initial as a state and as a named vector
  expect_equal(timeRelaxations(mc, x, initial = "a", timePoints = 0:3),
               timeRelaxations(mc, x, initial = c(b = 0, c = 0, a = 1),
                               timePoints = 0:3))
})

test_that("limits, order of time points and long lags", {
  f <- c(3, 2, 1)
  pi <- as.numeric(steadyStates(mc))
  expect_equal(unname(timeCorrelations(mc, x, timePoints = 2000)),
               sum(pi * f)^2, tolerance = 1e-10)
  expect_equal(unname(timeRelaxations(mc, x, timePoints = 1e6)), sum(pi * f),
               tolerance = 1e-10)
  # unsorted and repeated time points, and a lag reached by squaring
  expect_equal(unname(timeCorrelations(mc, x, timePoints = c(10, 0, 10, 100))),
               unname(timeCorrelations(mc, x, timePoints = c(10, 0, 10, 100) + 0)))
  direct <- function(t) { A <- diag(3); for (i in seq_len(t)) A <- A %*% P; sum(f * pi * (A %*% f)) }
  expect_equal(unname(timeCorrelations(mc, x, timePoints = c(100, 67))),
               c(direct(100), direct(67)), tolerance = 1e-12)
})

test_that("byrow = FALSE gives the same values", {
  mcc <- new("markovchain", transitionMatrix = t(P), byrow = FALSE)
  expect_equal(timeCorrelations(mcc, x, y, timePoints = tp),
               timeCorrelations(mc, x, y, timePoints = tp))
  expect_equal(timeRelaxations(mcc, x, timePoints = tp),
               timeRelaxations(mc, x, timePoints = tp))
})

test_that("timeCorrelations needs a unique stationary distribution", {
  red <- new("markovchain", transitionMatrix = diag(2),
             states = c("a", "b"))
  expect_error(timeCorrelations(red, c("a", "b")), "unique stationary")
  # relaxations are defined anyway
  expect_equal(unname(timeRelaxations(red, c("a", "a", "b"), timePoints = 0:2)),
               rep(1.5, 3))
})

test_that("arguments are validated", {
  expect_error(timeCorrelations(mc, c("a", "z")), "not in the chain")
  expect_error(timeCorrelations(mc, c("a", NA)), "missing")
  expect_error(timeCorrelations(mc, x, timePoints = -1), "timePoints")
  expect_error(timeCorrelations(mc, x, timePoints = 1.5), "timePoints")
  expect_error(timeRelaxations(mc, x, initial = c(0.5, 0.5)), "initial")
  expect_equal(timeCorrelations(mc, factor(x)), timeCorrelations(mc, x))
})
