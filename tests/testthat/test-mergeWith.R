library(testthat)
library(markovchain)

states2 <- c("a", "b")
mc1 <- new("markovchain", states = states2,
           transitionMatrix = matrix(c(0.9, 0.1, 0.1, 0.9), byrow = TRUE, nrow = 2,
                                     dimnames = list(states2, states2)))
mc2 <- new("markovchain", states = states2,
           transitionMatrix = matrix(c(0.5, 0.5, 0.5, 0.5), byrow = TRUE, nrow = 2,
                                     dimnames = list(states2, states2)))

## ---- basic convex combination ----------------------------------------------

test_that("mergeWith computes the convex combination entry-by-entry", {
  merged <- mergeWith(mc1, mc2, gamma = 0.25)
  expected <- 0.75 * mc1@transitionMatrix + 0.25 * mc2@transitionMatrix
  expect_equal(unclass(merged@transitionMatrix), unclass(expected),
               ignore_attr = TRUE)
  expect_true(merged@byrow)
})

test_that("gamma = 0 and gamma = 1 return the two inputs unchanged", {
  expect_equal(unclass(mergeWith(mc1, mc2, gamma = 0)@transitionMatrix),
               unclass(mc1@transitionMatrix), ignore_attr = TRUE)
  expect_equal(unclass(mergeWith(mc1, mc2, gamma = 1)@transitionMatrix),
               unclass(mc2@transitionMatrix), ignore_attr = TRUE)
})

test_that("the merged transition matrix is row-stochastic", {
  merged <- mergeWith(mc1, mc2, gamma = 0.5)
  expect_equal(as.numeric(rowSums(merged@transitionMatrix)), c(1, 1), tolerance = 1e-10)
})

## ---- states matched by name, not by position -------------------------------

test_that("mergeWith realigns states by name when other lists them in a different order", {
  # mc2_reordered describes exactly the same chain as mc2, but with the two
  # states listed in the opposite order and a different byrow convention.
  mc2_reordered <- new("markovchain", states = rev(states2), byrow = FALSE,
                       transitionMatrix = t(mc2@transitionMatrix[rev(states2), rev(states2)]))
  expect_equal(unclass(mc2_reordered@transitionMatrix[states2, states2]),
               unclass(mc2@transitionMatrix), ignore_attr = TRUE)

  merged1 <- mergeWith(mc1, mc2, gamma = 0.3)
  merged2 <- mergeWith(mc1, mc2_reordered, gamma = 0.3)

  expect_equal(unclass(merged1@transitionMatrix[states2, states2]),
               unclass(merged2@transitionMatrix[states2, states2]),
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("mergeWith rejects chains on different state sets", {
  other <- new("markovchain", states = c("a", "z"),
              transitionMatrix = matrix(c(0.5, 0.5, 0.5, 0.5), byrow = TRUE, nrow = 2,
                                        dimnames = list(c("a", "z"), c("a", "z"))))
  expect_error(mergeWith(mc1, other), "same set of state names")
})

## ---- validation --------------------------------------------------------------

test_that("mergeWith validates gamma", {
  expect_error(mergeWith(mc1, mc2, gamma = -0.1), "gamma")
  expect_error(mergeWith(mc1, mc2, gamma = 1.1), "gamma")
  expect_error(mergeWith(mc1, mc2, gamma = c(0.5, 0.5)), "gamma")
})

test_that("mergeWith(x, x, gamma) leaves a chain unchanged for any gamma", {
  merged <- mergeWith(mc1, mc1, gamma = 0.7)
  expect_equal(unclass(merged@transitionMatrix), unclass(mc1@transitionMatrix),
               tolerance = 1e-10, ignore_attr = TRUE)
})
