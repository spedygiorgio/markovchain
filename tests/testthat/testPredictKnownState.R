context("predict() and conditionalDistribution() validate the starting state (#145)")

test_that("an unknown or malformed starting state gives an informative error", {
  mc <- new("markovchain", states = c("a", "b", "c"),
            transitionMatrix = matrix(c(0.5, 0.5, 0,
                                        0.2, 0.8, 0,
                                        0.3, 0.3, 0.4), nrow = 3, byrow = TRUE))
  expect_error(predict(mc, newdata = "z"), "'z' is not one of the states")
  expect_error(predict(mc, newdata = c("a", "z"), n.ahead = 2), "newdata")
  expect_error(conditionalDistribution(mc, "z"), "'z' is not one of the states")
  expect_error(predict(mc, newdata = character(0)), "at least one state")
  expect_error(predict(mc, newdata = "a", n.ahead = 0), "positive integer")
  expect_error(predict(mc, newdata = "a", n.ahead = 1.5), "positive integer")
  expect_error(conditionalDistribution(mc, NA_character_), "must be a single state")
  expect_error(conditionalDistribution(mc, c("a", "b")), "must be a single state")
})

test_that("valid calls are unchanged", {
  mc <- new("markovchain", states = c("a", "b"),
            transitionMatrix = matrix(c(0.1, 0.9, 0.8, 0.2), nrow = 2, byrow = TRUE))
  expect_equal(predict(mc, newdata = "a", n.ahead = 3), c("b", "a", "b"))
  expect_equal(predict(mc, newdata = c("b", "a"), n.ahead = 1), "b")
  expect_equal(conditionalDistribution(mc, "a"), c(a = 0.1, b = 0.9))
  # a state seen only as the last element of a sequence is a valid state
  fit <- markovchainFit(c("a", "b", "a", "b", "c"))$estimate
  expect_length(predict(fit, newdata = "c", n.ahead = 2), 2)
})
