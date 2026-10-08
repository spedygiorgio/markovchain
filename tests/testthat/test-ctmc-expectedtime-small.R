# ExpectedTime() with two states: the matrix of the states other than j is
# 1 x 1 and must not lose its dimensions
test_that("ExpectedTime works on a CTMC with two states", {
  states <- c("sigma", "sigma_star")
  gen <- matrix(c(-3, 3,
                   1, -1), nrow = 2, byrow = TRUE,
                dimnames = list(states, states))
  mc <- new("ctmc", states = states, byrow = TRUE, generator = gen,
            name = "two states")
  for (useRCpp in c(TRUE, FALSE)) {
    expect_equal(ExpectedTime(mc, 1, 2, useRCpp = useRCpp), 1 / 3)
    expect_equal(ExpectedTime(mc, 2, 1, useRCpp = useRCpp), 1)
    expect_equal(ExpectedTime(mc, 1, 1, useRCpp = useRCpp), 0)
  }
  # generator stored by columns
  mcCol <- new("ctmc", states = states, byrow = FALSE, generator = t(gen),
               name = "two states by column")
  expect_equal(ExpectedTime(mcCol, 1, 2), 1 / 3)
})
