test_that("first-passage functions reject invalid dimensions and states", {
  states <- c("a", "b")
  mc <- new("markovchain", states = states,
    transitionMatrix = matrix(c(.5, .5, .4, .6), 2, byrow = TRUE,
      dimnames = list(states, states)))

  expect_error(firstPassage(mc, "missing", 2), "state")
  expect_error(firstPassage(mc, "a", 0), "positive integer")
  expect_error(firstPassageMultiple(mc, "a", character(), 2), "set")
  expect_error(firstPassageMultiple(mc, "a", "missing", 2), "set")
})

test_that("reward functions validate indices and vector lengths", {
  states <- c("a", "b")
  mc <- new("markovchain", states = states,
    transitionMatrix = matrix(c(.5, .5, .4, .6), 2, byrow = TRUE,
      dimnames = list(states, states)))

  expect_error(expectedRewards(mc, 1, 1), "one finite numeric value")
  expect_error(expectedRewards(mc, -1, c(1, 2)), "non-negative integer")
  expect_error(expectedRewardsBeforeHittingA(
    mc, A = "b", state = "missing", rewards = c(1, 2), n = 2), "state")
  expect_error(expectedRewardsBeforeHittingA(
    mc, A = "b", state = "b", rewards = c(1, 2), n = 2), "must not belong")
})

test_that("higher-order helpers reject unsafe orders", {
  expect_error(seq2freqProb(character()), "must not be empty")
  expect_error(seq2matHigh(c("a", "b"), -1), "order")
  expect_error(seq2matHigh(c("a", "b"), 2), "order")
  expect_error(fitHigherOrder(c("a", "b"), 2), "order")
})

test_that("visit and simulation functions validate sizes", {
  states <- c("a", "b")
  mc <- new("markovchain", states = states,
    transitionMatrix = matrix(c(.5, .5, .4, .6), 2, byrow = TRUE,
      dimnames = list(states, states)))

  expect_error(noofVisitsDist(mc, 0, "a"), "positive integer")
  expect_error(noofVisitsDist(mc, 2, "missing"), "initial state")
  expect_error(rmarkovchain(0, mc), "positive integer")
})

test_that("committorAB validates complete sets", {
  states <- c("a", "b", "c")
  mc <- new("markovchain", states = states,
    transitionMatrix = diag(3),
    name = "identity")

  expect_error(committorAB(mc, c(0, 2), 3), "set A")
  expect_error(committorAB(mc, c(1, 2), c(2, 3)), "disjoint")
})

## ---- hardening round (1.2.1): arguments that used to be accepted silently ----

test_that("multinomialConfidenceIntervals validates the confidence level", {
  P <- matrix(c(0.5, 0.5, 0.25, 0.75), 2, byrow = TRUE, dimnames = list(c("a", "b"), c("a", "b")))
  N <- matrix(c(5, 5, 3, 9), 2, byrow = TRUE, dimnames = dimnames(P))
  for (bad in list(0, 1, 2, -0.1, NA_real_, Inf, c(0.9, 0.95), "0.95")) {
    expect_error(multinomialConfidenceIntervals(P, N, confidencelevel = bad), "confidencelevel")
  }
  ci <- multinomialConfidenceIntervals(P, N, confidencelevel = 0.9)
  expect_true(all(ci$lowerEndpointMatrix <= ci$upperEndpointMatrix))
})

test_that("priorDistribution and inferHyperparam reject NA hyperparameters and bad scales", {
  st <- c("a", "b")
  P <- matrix(0.5, 2, 2, dimnames = list(st, st))
  H <- matrix(2, 2, 2, dimnames = list(st, st))
  H[1, 2] <- NA
  expect_error(priorDistribution(P, H), "greater than or equal to 1")
  expect_error(inferHyperparam(P, scale = c(10, -5)), "positive")
  expect_error(inferHyperparam(P, scale = c(10, NA)), "positive")
  expect_error(inferHyperparam(P, scale = c(10, 0)), "positive")
  expect_silent(inferHyperparam(P, scale = c(10, 10)))
})

test_that("absurdly large n gives a clear error", {
  mc <- new("markovchain", states = c("a", "b"),
            transitionMatrix = matrix(c(0.5, 0.5, 0.5, 0.5), 2, byrow = TRUE))
  expect_error(firstPassage(mc, "a", 2^31), "too large")
  expect_error(markovchainSequence(2^31, mc, t0 = "a"), "must not exceed")
})
