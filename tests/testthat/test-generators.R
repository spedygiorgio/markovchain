library(testthat)
library(markovchain)

## ---- birthDeath -------------------------------------------------------------

test_that("birthDeath builds the expected tridiagonal row-stochastic matrix", {
  bd <- birthDeath(p = c(0.3, 0.4, 0.5), q = c(0.2, 0.3, 0.1))
  expected <- matrix(c(
    0.7, 0.3, 0.0, 0.0,
    0.2, 0.4, 0.4, 0.0,
    0.0, 0.3, 0.2, 0.5,
    0.0, 0.0, 0.1, 0.9
  ), nrow = 4, byrow = TRUE, dimnames = list(as.character(1:4), as.character(1:4)))
  expect_equal(unclass(bd@transitionMatrix), unclass(expected), ignore_attr = TRUE)
  expect_equal(as.numeric(rowSums(bd@transitionMatrix)), rep(1, 4))
  expect_true(bd@byrow)
})

test_that("birthDeath accepts custom state names", {
  bd <- birthDeath(p = c(0.5), q = c(0.5), states = c("low", "high"))
  expect_equal(states(bd), c("low", "high"))
})

test_that("birthDeath validates its arguments", {
  expect_error(birthDeath(p = c(1.5), q = c(0.5)), "p must be")
  expect_error(birthDeath(p = c(0.5, 0.5), q = c(0.5)), "q must be")
  # 3-state chain: interior state 2 has p[2]+q[1] = 0.5+0.6 = 1.1 > 1.
  expect_error(birthDeath(p = c(0.6, 0.5), q = c(0.6, 0.5)), "diagonal")
  expect_error(birthDeath(p = c(0.5), q = c(0.5), states = c("a")), "states")
})

## ---- gamblersRuin -------------------------------------------------------------

test_that("gamblersRuin has absorbing ends and the expected interior structure", {
  ruin <- gamblersRuin(upperBound = 5, prob = 0.4)
  expect_equal(states(ruin), as.character(0:5))
  expect_equal(sort(absorbingStates(ruin)), c("0", "5"))
  expect_equal(as.numeric(rowSums(ruin@transitionMatrix)), rep(1, 6))
  expect_equal(ruin@transitionMatrix["2", "1"], 0.6)
  expect_equal(ruin@transitionMatrix["2", "3"], 0.4)
})

test_that("gamblersRuin matches the classical ruin probability formula", {
  upperBound <- 6
  prob <- 0.35
  ruin <- gamblersRuin(upperBound = upperBound, prob = prob)
  ap <- absorptionProbabilities(ruin)
  r <- (1 - prob) / prob
  i <- 1:(upperBound - 1)
  theoretical <- (r^i - r^upperBound) / (1 - r^upperBound)
  # absorptionProbabilities returns probability of ending in each absorbing
  # state; column "0" is the ruin probability.
  expect_equal(as.numeric(ap[as.character(i), "0"]), theoretical, tolerance = 1e-6)
})

test_that("gamblersRuin validates its arguments", {
  expect_error(gamblersRuin(upperBound = 1, prob = 0.5), "upperBound")
  expect_error(gamblersRuin(upperBound = 5, prob = 1.5), "prob")
})

## ---- urnModel -------------------------------------------------------------

test_that("urnModel (Ehrenfest) has Binomial(balls, 0.5) stationary distribution", {
  ehr <- urnModel(balls = 6)
  expect_equal(as.numeric(rowSums(ehr@transitionMatrix)), rep(1, 7))
  pi <- as.numeric(steadyStates(ehr))
  expect_equal(pi, dbinom(0:6, 6, 0.5), tolerance = 1e-8)
})

test_that("urnModel is reflecting (not absorbing) at the boundaries", {
  ehr <- urnModel(balls = 4)
  expect_equal(ehr@transitionMatrix["0", "1"], 1)
  expect_equal(ehr@transitionMatrix["4", "3"], 1)
  expect_length(absorbingStates(ehr), 0L)
})

test_that("urnModel validates its arguments", {
  expect_error(urnModel(balls = 0), "balls")
  expect_error(urnModel(balls = -1), "balls")
})

## ---- toBoundedChain -------------------------------------------------------------

bd <- birthDeath(p = c(0.3, 0.4, 0.5), q = c(0.2, 0.3, 0.1))

test_that("toBoundedChain('absorbing') makes the first and last state absorbing", {
  absorbed <- toBoundedChain(bd, "absorbing")
  expect_equal(sort(absorbingStates(absorbed)), c("1", "4"))
  # Interior rows are untouched.
  expect_equal(absorbed@transitionMatrix[2:3, ], bd@transitionMatrix[2:3, ])
})

test_that("toBoundedChain('reflecting') forces deterministic bounce-back", {
  reflected <- toBoundedChain(bd, "reflecting")
  expect_equal(as.numeric(reflected@transitionMatrix[1, ]), c(0, 1, 0, 0))
  expect_equal(as.numeric(reflected@transitionMatrix[4, ]), c(0, 0, 1, 0))
})

test_that("toBoundedChain(numeric) interpolates between absorbing and reflecting", {
  expect_equal(unclass(toBoundedChain(bd, 0)@transitionMatrix),
               unclass(toBoundedChain(bd, "absorbing")@transitionMatrix),
               ignore_attr = TRUE)
  expect_equal(unclass(toBoundedChain(bd, 1)@transitionMatrix),
               unclass(toBoundedChain(bd, "reflecting")@transitionMatrix),
               ignore_attr = TRUE)

  semi <- toBoundedChain(bd, 0.25)
  expect_equal(as.numeric(semi@transitionMatrix[1, ]), c(0.75, 0.25, 0, 0))
  expect_equal(as.numeric(rowSums(semi@transitionMatrix)), rep(1, 4))
})

test_that("toBoundedChain validates its arguments", {
  expect_error(toBoundedChain(bd, "invalid"), "boundaryCondition")
  expect_error(toBoundedChain(bd, -0.1), "boundaryCondition")
  expect_error(toBoundedChain(bd, 1.1), "boundaryCondition")

  one_state <- new("markovchain", states = "a",
                   transitionMatrix = matrix(1, 1, 1, dimnames = list("a", "a")))
  expect_error(toBoundedChain(one_state, "absorbing"), "at least 2 states")
})

## ---- tauchen and rouwenhorst -------------------------------------------------------------

test_that("tauchen and rouwenhorst return row-stochastic chains with matching grids", {
  out_t <- tauchen(alpha = 0, sigma = 1, rho = 0.5, size = 7)
  out_r <- rouwenhorst(alpha = 0, sigma = 1, rho = 0.5, size = 7)

  expect_equal(as.numeric(rowSums(out_t$chain@transitionMatrix)), rep(1, 7), tolerance = 1e-10)
  expect_equal(as.numeric(rowSums(out_r$chain@transitionMatrix)), rep(1, 7), tolerance = 1e-10)
  expect_length(out_t$states, 7)
  expect_length(out_r$states, 7)
  expect_true(all(diff(out_t$states) > 0))
  expect_true(all(diff(out_r$states) > 0))
})

test_that("rouwenhorst exactly reproduces the AR(1) unconditional variance and lag-1 autocorrelation", {
  alpha <- 1; sigma <- 2; rho <- 0.9
  out <- rouwenhorst(alpha = alpha, sigma = sigma, rho = rho, size = 9)
  pi <- as.numeric(steadyStates(out$chain))
  y <- out$states
  mu <- sum(pi * y)
  varY <- sum(pi * (y - mu)^2)
  theoreticalVar <- sigma^2 / (1 - rho^2)
  expect_equal(mu, alpha, tolerance = 1e-8)
  expect_equal(varY, theoreticalVar, tolerance = 1e-8)

  condMean <- as.numeric(out$chain@transitionMatrix %*% y)
  empiricalRho <- sum(pi * (y - mu) * (condMean - mu)) / varY
  expect_equal(empiricalRho, rho, tolerance = 1e-8)
})

test_that("rouwenhorst stays accurate for highly persistent (near unit-root) processes where tauchen drifts", {
  alpha <- 0; sigma <- 1; rho <- 0.98
  theoreticalVar <- sigma^2 / (1 - rho^2)

  out_r <- rouwenhorst(alpha = alpha, sigma = sigma, rho = rho, size = 15)
  pi_r <- as.numeric(steadyStates(out_r$chain))
  var_r <- sum(pi_r * (out_r$states - sum(pi_r * out_r$states))^2)
  expect_equal(var_r, theoreticalVar, tolerance = 1e-6)

  out_t <- tauchen(alpha = alpha, sigma = sigma, rho = rho, size = 15)
  pi_t <- as.numeric(steadyStates(out_t$chain))
  var_t <- sum(pi_t * (out_t$states - sum(pi_t * out_t$states))^2)
  # tauchen's fixed k*sigma_y grid width is not designed to match the
  # unconditional variance exactly; rouwenhorst's is markedly closer here.
  expect_gt(abs(var_t - theoreticalVar), abs(var_r - theoreticalVar))
})

test_that("tauchen and rouwenhorst validate their arguments", {
  expect_error(tauchen(alpha = 0, sigma = -1, rho = 0.5, size = 5), "sigma")
  expect_error(tauchen(alpha = 0, sigma = 1, rho = 1, size = 5), "rho")
  expect_error(tauchen(alpha = 0, sigma = 1, rho = 0.5, size = 1), "size")
  expect_error(tauchen(alpha = 0, sigma = 1, rho = 0.5, size = 5, k = -1), "k")

  expect_error(rouwenhorst(alpha = 0, sigma = 0, rho = 0.5, size = 5), "sigma")
  expect_error(rouwenhorst(alpha = 0, sigma = 1, rho = -1, size = 5), "rho")
  expect_error(rouwenhorst(alpha = 0, sigma = 1, rho = 0.5, size = 1), "size")
})
