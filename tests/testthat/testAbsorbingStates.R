test_that("markovchainFit supports explicit absorbing states", {
  c1 <- c("a", "b", "c", "c", "e")
  c2 <- c("a", "b", "d", "e")
  c3 <- c("a", "c", "b", "c", "d")
  c4 <- c("a", "b", "b", "d", "b", "c", "d", "e")
  c5 <- c("a", "c", "c", "d", "d")
  c6 <- c("a", "c", "d", "d", "b", "b", "e")
  journeys <- list(c1, c2, c3, c4, c5, c6)

  fit <- markovchainFit(journeys, absorbingStates = "e")
  p <- fit$estimate@transitionMatrix

  expect_equal(p["e", ], c(a = 0, b = 0, c = 0, d = 0, e = 1))
  expect_equal(sum(p["e", ]), 1)
  expect_equal(p["a", "b"], 0.5)
  expect_equal(p["d", "e"], 1 / 3)
})

test_that("multiple absorbing states are supported", {
  journeys <- list(
    c("a", "b", "end"),
    c("a", "c", "death")
  )

  fit <- markovchainFit(
    journeys,
    absorbingStates = c("end", "death")
  )
  p <- fit$estimate@transitionMatrix

  expect_equal(p["end", ], c(a = 0, b = 0, death = 0, end = 1))
  expect_equal(p["death", ], c(a = 0, b = 0, death = 1, end = 0))
})

test_that("absorbing states respect byrow = FALSE", {
  journeys <- list(
    c("a", "b", "end"),
    c("a", "b", "end"),
    c("a", "c", "end")
  )

  fit <- markovchainFit(
    journeys,
    byrow = FALSE,
    absorbingStates = "end"
  )
  p <- fit$estimate@transitionMatrix

  expect_equal(p[, "end"], c(a = 0, b = 0, c = 0, end = 1))
  expect_equal(sum(p[, "end"]), 1)
  expect_equal(p["b", "a"], 1)
})

test_that("declared absorbing states cannot have observed outgoing transitions", {
  journeys <- list(
    c("a", "end", "b"),
    c("a", "b")
  )

  expect_error(
    markovchainFit(journeys, absorbingStates = "end"),
    "observed outgoing transitions"
  )
})

test_that("absorbingStates is only available for MLE", {
  expect_error(
    markovchainFit(c("a", "b", "a"), method = "laplace", absorbingStates = "b"),
    "only supported with method = \"mle\""
  )
})

test_that("empty absorbingStates preserves the existing fitting path", {
  sequence <- c("a", "b", "a", "a", "b")

  fit_default <- markovchainFit(sequence)
  fit_empty <- markovchainFit(sequence, absorbingStates = character())

  expect_equal(
    fit_empty$estimate@transitionMatrix,
    fit_default$estimate@transitionMatrix
  )
})
