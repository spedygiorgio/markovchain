library(markovchain)

#create basic markov chains
require(matlab)
mathematicaMatr <- zeros(5)
mathematicaMatr[1,] <- c(0, 1/3, 0, 2/3, 0)
mathematicaMatr[2,] <- c(1/2, 0, 0, 0, 1/2)
mathematicaMatr[3,] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[4,] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[5,] <- c(0, 0, 0, 0, 1)
statesNames <- letters[1:5]
mathematicaMc <- new("markovchain", transitionMatrix = mathematicaMatr,
                     name = "Mathematica MC", states = statesNames)

# ####end creating DTMC
# test .gcdRcpp function, .commclassesKernelRcpp function, .commStatesFinderRcpp function
context("Optimization of functions")

test_that("Optimized functions should satisfy", 
          {
            expect_equal(.gcdRcpp(9, 12), 3) # .gcdRcpp function is also tested in testPeriod.R
            expect_equal(.commStatesFinderRcpp(mathematicaMatr), matrix(c(1,1,1,1,1,
                                                                          1,1,1,1,1,
                                                                          0,0,1,1,0,
                                                                          0,0,1,1,0,
                                                                          0,0,0,0,1
                                                                          ), nrow=5, byrow=T
            ))
          })

