library(markovchain)

#create basic markov chains
markov1<-new("markovchain", states=c("a","b","c"), transitionMatrix=
               matrix(c(0.2,0.5,0.3,
                        0,1,0,
                        0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(c("a","b","c"),
                                                                       c("a","b","c"))
               ))

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
####end creating DTMC
context("Basic DTMC proprieties")

test_that("States are those that should be", {
  expect_equal(absorbingStates(markov1), "b")
  expect_equal(transientStates(markov1), c("a","c"))
  expect_equal(is.irreducible(mathematicaMc),FALSE)
  expect_equal(transientStates(mathematicaMc), c("a","b"))
  expect_equal(is.accessible(mathematicaMc, "a", "c"),TRUE)
  expect_equal(.canonicForm(mathematicaMc)@transitionMatrix, .canonicFormRcpp(mathematicaMc)@transitionMatrix) 
  expect_equal(summary(mathematicaMc), list(closedClasses = list(c("c", "d"), c("e")), 
                                            transientClasses = list(c("a", "b"))))
})

###testing proper conversion of objects
context("Conversion of objects")
provaMatr2Mc<-as(mathematicaMatr,"markovchain")

test_that("Conversion of objects", 
          {
            expect_equal(class(provaMatr2Mc)=="markovchain",TRUE)
          })

###perform some fitting
sequence1<-c("a", "b", "a", "a", "a")
sequence2<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
mcFit<-markovchainFit(data=sequence1,byrow=FALSE)
test_that("Fit should satisfy", {
  expect_equal((mcFit["logLikelihood"])[[1]], log(1/3) + 2*log(2/3))
  expect_equal(markovchainFit(data=sequence2, method="bootstrap")["confidenceInterval"]
               [[1]]["confidenceLevel"][[1]], 0.95)
})

data(rain)
sequs<-rain$rain
# mcBoot<-markovchainFit(data = sequs,nboot = 10,method="bootstrap") #ok
# mcBoot2<-markovchainFit(data = sequs,nboot = 200,method="bootstrap") #ok but slower
mcBoot<-markovchainFit(data = sequs[1:100],nboot = 100,method="bootstrap",parallel=TRUE) # ok
# mcBoot<-markovchainFit(data = sequs,nboot = 100,method="bootstrap",parallel=TRUE)

### MAP fit function tests
data1 <- c("a", "b", "a", "c", "a", "b", "a", "b", "c", "b", "b", "a", "b")
data2 <- c("c", "a", "b")

test_that("MAP fits must satisfy", {
  expect_identical(markovchainFit(data1, method = "map")$estimate@transitionMatrix, 
                   markovchainFit(data1, method = "mle")$estimate@transitionMatrix)
  
  expect_identical(markovchainFit(data1, method = "map")$estimate@transitionMatrix, 
                   matrix(c(0.0, 0.6, 0.5, 
                           0.8, 0.2, 0.5, 
                           0.2, 0.2, 0.0), nrow = 3, 
                          dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))
  
  expect_identical(markovchainFit(data1, method = "map", hyperparam = 
                                    matrix(c(2, 1, 3, 
                                             4, 5, 2, 
                                             2, 2, 1), nrow = 3, 
                                           dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))$estimate@transitionMatrix, 
                   matrix(c(1/10, 3/10, 3/5, 
                            7/10, 5/10, 2/5, 
                            2/10, 2/10, 0), nrow = 3, 
                          dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))
})

test_that("predictiveDistribution must satisfy", {
  expect_equal(predictiveDistribution(data1, character()), 0)
  
  expect_equal(predictiveDistribution(data1, data2, hyperparam = 
                                            matrix(c(2, 1, 3, 
                                                     4, 5, 2, 
                                                     2, 2, 1), nrow = 3, 
                                                   dimnames = list(c("a", "b", "c"), c("a", "b", "c")))), 
                   log(4 / 13))
})

test_that("inferHyperparam must satisfy", {
  expect_identical(inferHyperparam(data = data1)$dataInference, 
                   matrix(c(1, 4, 2, 
                            5, 2, 2, 
                            2, 2, 1), nrow = 3, 
                          dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))
  
  expect_identical(inferHyperparam(transMatr = 
                                     matrix(c(0.0, 0.6, 0.5, 
                                              0.8, 0.2, 0.5, 
                                              0.2, 0.2, 0.0), nrow = 3, 
                                            dimnames = list(c("a", "b", "c"), c("a", "b", "c"))),
                                   scale = c(10, 10, 10))$scaledInference, 
                   matrix(c(0, 6, 5, 
                            8, 2, 5, 
                            2, 2, 0), nrow = 3, 
                          dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))
})
