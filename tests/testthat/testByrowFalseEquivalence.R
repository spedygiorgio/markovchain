context("Column-stochastic storage (byrow = FALSE) gives the same results as row storage (#148)")

.byrowFixtures <- function() {
  states <- c("a", "b", "c")
  P <- matrix(c(0.5, 0.3, 0.2,
                0.1, 0.6, 0.3,
                0.2, 0.2, 0.6),
              nrow = 3, byrow = TRUE, dimnames = list(states, states))
  list(
    P = P,
    byRow = new("markovchain", transitionMatrix = P, byrow = TRUE),
    byCol = new("markovchain", transitionMatrix = t(P), byrow = FALSE)
  )
}

test_that("sequence generation follows the outgoing probabilities for both orientations", {
  fx <- .byrowFixtures()
  for (useRCpp in c(TRUE, FALSE)) {
    set.seed(11)
    sRow <- markovchainSequence(2000, fx$byRow, t0 = "a", useRCpp = useRCpp)
    set.seed(11)
    sCol <- markovchainSequence(2000, fx$byCol, t0 = "a", useRCpp = useRCpp)
    expect_identical(sRow, sCol)
  }
  set.seed(12)
  big <- markovchainSequence(40000, fx$byCol, t0 = "a", include.t0 = TRUE)
  emp <- prop.table(table(factor(head(big, -1), levels = rownames(fx$P)),
                          factor(big[-1], levels = colnames(fx$P))), 1)
  expect_equal(as.numeric(emp), as.numeric(fx$P), tolerance = 0.03)
})

test_that("rmarkovchain gives identical draws for both orientations", {
  fx <- .byrowFixtures()
  set.seed(3); xRow <- rmarkovchain(200, fx$byRow, t0 = "a")
  set.seed(3); xCol <- rmarkovchain(200, fx$byCol, t0 = "a")
  expect_identical(xRow, xCol)
  list_row <- new("markovchainList", markovchains = list(fx$byRow, fx$byRow))
  list_col <- new("markovchainList", markovchains = list(fx$byCol, fx$byCol))
  set.seed(5); yRow <- rmarkovchain(20, list_row, t0 = "a")
  set.seed(5); yCol <- rmarkovchain(20, list_col, t0 = "a")
  expect_identical(yRow, yCol)
})

test_that("committorAB, expectedRewards and related quantities do not depend on the orientation", {
  fx <- .byrowFixtures()
  expect_equal(committorAB(fx$byRow, 1, 3), committorAB(fx$byCol, 1, 3))
  # hand computation: q(b) = 0.1 + 0.6 q(b)  =>  q(b) = 0.25
  expect_equal(unname(committorAB(fx$byCol, 1, 3))[2], 0.25)
  expect_equal(expectedRewards(fx$byRow, 2, c(1, 2, 3)),
               expectedRewards(fx$byCol, 2, c(1, 2, 3)))
  expect_equal(expectedRewardsBeforeHittingA(fx$byRow, "c", "a", c(1, 2, 3), 3),
               expectedRewardsBeforeHittingA(fx$byCol, "c", "a", c(1, 2, 3), 3))
  expect_equal(firstPassage(fx$byRow, "a", 4), firstPassage(fx$byCol, "a", 4))
  expect_equal(firstPassageMultiple(fx$byRow, "a", c("b", "c"), 4),
               firstPassageMultiple(fx$byCol, "a", c("b", "c"), 4))
  expect_equal(noofVisitsDist(fx$byRow, 5, "a"), noofVisitsDist(fx$byCol, 5, "a"))
  expect_equal(is.stochasticallyMonotone(fx$byRow), is.stochasticallyMonotone(fx$byCol))
})

test_that("lump and autoLump work for column-stochastic chains", {
  states <- letters[1:5]
  B <- matrix(c(.9, .1, 0, 0, 0,
                .1, .9, 0, 0, 0,
                0, 0, .5, .3, .2,
                0, 0, .3, .4, .3,
                0, 0, .2, .3, .5),
              nrow = 5, byrow = TRUE, dimnames = list(states, states))
  byRow <- new("markovchain", transitionMatrix = B, byrow = TRUE)
  byCol <- new("markovchain", transitionMatrix = t(B), byrow = FALSE)
  partition <- list(X = c("a", "b"), Y = c("c", "d", "e"))
  expect_equal(lump(byRow, partition)@transitionMatrix,
               lump(byCol, partition)@transitionMatrix)
  expect_equal(autoLump(byRow, 2)$partition, autoLump(byCol, 2)$partition)
})
