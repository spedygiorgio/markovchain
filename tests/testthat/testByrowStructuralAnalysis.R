context("Structural analysis with column-stochastic matrices")

test_that("column-stochastic matrices preserve recurrent and transient classes", {
  states <- as.character(0:4)
  P <- matrix(
    c(
      1, 0, 0, 0, 0,
      0.5, 0.5, 0, 0, 0,
      0, 0.5, 0, 0.5, 0,
      0, 0, 0.5, 0.5, 0,
      0, 0, 0, 0, 1
    ),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(states, states)
  )

  mc_byrow <- new(
    "markovchain",
    states = states,
    transitionMatrix = P,
    byrow = TRUE
  )

  mc_bycol <- new(
    "markovchain",
    states = states,
    transitionMatrix = t(P),
    byrow = FALSE
  )

  expect_equal(recurrentClasses(mc_byrow), recurrentClasses(mc_bycol))
  expect_equal(transientClasses(mc_byrow), transientClasses(mc_bycol))
  expect_equal(recurrentStates(mc_byrow), recurrentStates(mc_bycol))
  expect_equal(transientStates(mc_byrow), transientStates(mc_bycol))

  expect_equal(recurrentClasses(mc_bycol), list("0", "4"))
  expect_equal(transientClasses(mc_bycol), list("1", c("2", "3")))
  expect_equal(recurrentStates(mc_bycol), c("0", "4"))
  expect_equal(transientStates(mc_bycol), c("1", "2", "3"))
})

test_that("column-stochastic matrices preserve reachability", {
  states <- c("a", "b", "c")
  P <- matrix(
    c(
      0, 1, 0,
      0, 0, 1,
      0, 0, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(states, states)
  )

  mc_byrow <- new(
    "markovchain",
    states = states,
    transitionMatrix = P,
    byrow = TRUE
  )

  mc_bycol <- new(
    "markovchain",
    states = states,
    transitionMatrix = t(P),
    byrow = FALSE
  )

  expect_equal(is.accessible(mc_byrow), is.accessible(mc_bycol))
  expect_true(is.accessible(mc_bycol)["a", "c"])
  expect_true(is.accessible(mc_bycol)["b", "c"])
  expect_false(is.accessible(mc_bycol)["c", "a"])
})
