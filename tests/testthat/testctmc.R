library(markovchain)

# Example from the book Markovchains, J. R. Norris, Cambridge University Press
states <- c("a","b","c","d")
byRow <- TRUE
gen <- matrix(data = c(-1, 1/2, 1/2, 0, 1/4, -1/2, 0, 1/4, 1/6, 0, -1/3, 1/6, 0, 0, 0, 0),
                  nrow = 4,byrow = byRow, dimnames = list(states,states))
ctmc <- new("ctmc",states = states, byrow = byRow, generator = gen, name = "testctmc")

test_that("Check Expected hitting time from one state to another",{
  expect_equal(ExpectedTime(ctmc,1,4),7)
  expect_equal(ExpectedTime(ctmc,2,4),5.5)
})