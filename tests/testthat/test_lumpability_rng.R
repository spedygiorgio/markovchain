context("Lumpability RNG handling")

test_that("autoLump does not alter the caller's RNG state", {
  states <- c("A", "B", "C", "D")
  P <- matrix(
    c(0.8, 0.1, 0.05, 0.05,
      0.1, 0.8, 0.05, 0.05,
      0.05, 0.05, 0.8, 0.1,
      0.05, 0.05, 0.1, 0.8),
    nrow = 4, byrow = TRUE,
    dimnames = list(states, states)
  )
  mc <- new("markovchain", states = states, transitionMatrix = P)

  set.seed(20260819)
  expected <- runif(5)

  set.seed(20260819)
  autoLump(mc, 2)
  actual <- runif(5)

  expect_equal(actual, expected)
})

test_that("autoLump preserves an absent RNG state", {
  states <- c("A", "B", "C", "D")
  P <- matrix(
    c(0.8, 0.1, 0.05, 0.05,
      0.1, 0.8, 0.05, 0.05,
      0.05, 0.05, 0.8, 0.1,
      0.05, 0.05, 0.1, 0.8),
    nrow = 4, byrow = TRUE,
    dimnames = list(states, states)
  )
  mc <- new("markovchain", states = states, transitionMatrix = P)

  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  autoLump(mc, 2)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})
