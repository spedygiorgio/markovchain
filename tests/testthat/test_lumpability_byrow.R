context("Lumpability with column-stochastic matrices")

test_that("lumpability is invariant to matrix orientation", {
  states <- c("rainy", "nice", "snowy")
  P <- matrix(
    c(0.5, 0.25, 0.25,
      0.5, 0.00, 0.50,
      0.25, 0.25, 0.50),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(states, states)
  )
  partition <- list(
    Bad_Weather = c("rainy", "snowy"),
    Nice_Weather = "nice"
  )

  mc_byrow <- new("markovchain", states = states,
                  transitionMatrix = P, byrow = TRUE)
  mc_bycol <- new("markovchain", states = states,
                  transitionMatrix = t(P), byrow = FALSE)

  expect_true(is.lumpable(mc_byrow, partition))
  expect_true(is.lumpable(mc_bycol, partition))

  lumped_byrow <- lump(mc_byrow, partition)
  lumped_bycol <- lump(mc_bycol, partition)

  expect_equal(lumped_bycol@transitionMatrix,
               lumped_byrow@transitionMatrix,
               tolerance = 1e-12)
  expect_true(lumped_bycol@byrow)
})

test_that("lumpability methods are exported S4 methods", {
  expect_true(isTRUE("is.lumpable" %in% getNamespaceExports("markovchain")))
  expect_true(isTRUE("lump" %in% getNamespaceExports("markovchain")))
  expect_true(isTRUE("autoLump" %in% getNamespaceExports("markovchain")))
})

test_that("partition names are validated", {
  mc <- new("markovchain", states = c("A", "B", "C"),
            transitionMatrix = matrix(1 / 3, 3, 3,
                                      dimnames = list(c("A", "B", "C"),
                                                      c("A", "B", "C"))))

  expect_error(is.lumpable(mc, list(c("A"), c("B", "C"))),
               "named list")
  expect_error(is.lumpable(mc, list(Group = c("A"), Group = c("B", "C"))),
               "unique")
})
