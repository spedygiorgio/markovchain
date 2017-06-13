library(markovchain)

context("Checking that ExpectedTime function works as expected")
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




context("Checking that probabilityatT function works as expected")
# TESTS for probabilityatT function
# Example taken from the book INTRODUCTION TO STOCHASTIC PROCESSES WITH R, ROBERT P. DOBROW, Wiley



states <- c("a","b","c","d","e")


# taken exactly from book
ansMatrix <- matrix(data = c(0.610, 0.290, 0.081, 0.016, 0.003,
                              0.232, 0.443, 0.238, 0.071, 0.017,
                              0.052, 0.190, 0.435, 0.238, 0.085,
                              0.008, 0.045, 0.191, 0.446, 0.310,
                              0.001, 0.008, 0.054, 0.248, 0.688),nrow = 5,byrow = T,dimnames = list(states,states))

byRow <- TRUE
gen <- matrix(c(-1/4,1/4,0,0,0,1/5,-9/20,1/4,0,0,0,1/5,-9/20,1/4,0,0,0,1/5,-9/20,1/4,0,0,0,1/5,-1/5),
              nrow=5,byrow=byRow, dimnames = list(states,states))

ctmc <- new("ctmc",states = states, byrow = byRow, generator = gen, name = "testctmc")

test_that("Check probabilityatT using a ctmc object:",{
  expect_equal(round(probabilityatT(ctmc,2.5),3),ansMatrix)
})


### Adds tests for impreciseprobabilityatT function
context("Checking that impreciseprobabilityatT function works as expected:")
states <- c("n","y")
Q <- matrix(c(-1,1,1,-1),nrow = 2,byrow = T,dimnames = list(states,states))
range <- matrix(c(1/52,3/52,1/2,2),nrow = 2,byrow = 2)
name <- "testictmc"
ictmc <- new("ictmc",states = states,Q = Q,range = range,name = name)


test_that("Check impreciseProbabilityatT function using an ictmc object:",{
  expect_equal(round(impreciseProbabilityatT(ictmc,2,0,1,error = 10^-3),4),c(0.0083,0.1410))
})
