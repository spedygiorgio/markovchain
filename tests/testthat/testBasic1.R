#library(markovchain)

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
  expect_equal(.recurrentClassesRcpp(mathematicaMc), list(c("c", "d"), c("e")))
  expect_equal(summary(mathematicaMc), list(closedClasses = list(c("c", "d"), c("e")), 
                                            recurrentClasses = list(c("c", "d"), c("e")),
                                            transientClasses = list(c("a", "b"))))
})

###testing proper conversion of objects
context("Conversion of objects")
provaMatr2Mc<-as(mathematicaMatr,"markovchain")

test_that("Conversion of objects", 
          {
            expect_equal(class(provaMatr2Mc)=="markovchain",TRUE)
          })

### Markovchain Fitting
sequence1 <- c("a", "b", "a", "a", "a")
sequence2 <- c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
mcFit <- markovchainFit(data = sequence1, byrow = FALSE, sanitize = TRUE)
mcFit2 <- markovchainFit(c("a","b","a","b"), sanitize = TRUE)

test_that("Fit should satisfy", {
  expect_equal((mcFit["logLikelihood"])[[1]], log(1/3) + 2*log(2/3))
  expect_equal(markovchainFit(data = sequence2, method = "bootstrap")["confidenceInterval"]
               [[1]]["confidenceLevel"][[1]], 0.95)
  expect_equal(mcFit2$confidenceInterval$upperEndpointMatrix, matrix(c(0,1,1,0), nrow = 2, byrow = TRUE,
                                                dimnames = list(c("a", "b"), c("a", "b"))))
})

### Markovchain Fitting for bigger markov chain
bigseq <- rep(c("a", "b", "c"), 500000)
bigmcFit <- markovchainFit(bigseq)

test_that("MC Fit for large sequence 1", {
  expect_equal(bigmcFit$logLikelihood, 0)
  expect_equal(bigmcFit$confidenceInterval$confidenceLevel, 0.95)
  expect_equal(bigmcFit$estimate@transitionMatrix, bigmcFit$confidenceInterval$upperEndpointMatrix)
})

bigmcFit <- markovchainFit(bigseq, sanitize = TRUE)

test_that("MC Fit for large sequence 2", {
  expect_equal(bigmcFit$logLikelihood, 0)
  expect_equal(bigmcFit$confidenceInterval$confidenceLevel, 0.95)
  expect_equal(bigmcFit$estimate@transitionMatrix, bigmcFit$confidenceInterval$upperEndpointMatrix)
})

### Markovchain Fitting For dataframe or matrix as an input
matseq <- matrix(c("a", "b", "c", "a","b", "c"), nrow = 2, byrow = T)

# for matrix as input

test_that("Markovchain Fit for matrix as input", {

# for matrix as input    
  
  expect_equal(markovchainFit(matseq)$estimate@transitionMatrix, 
               matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), nrow = 3, 
                      byrow = TRUE, dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))
  
  expect_equal(markovchainFit(matseq, sanitize = TRUE)$estimate@transitionMatrix, 
               matrix(c(0, 1, 0, 0, 0, 1, 1/3, 1/3, 1/3), nrow = 3, 
                      byrow = TRUE, dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))
  
# for data frame as input
    expect_equal(markovchainFit(as.data.frame(matseq))$estimate@transitionMatrix, 
                 matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), nrow = 3, 
                        byrow = TRUE, dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))
    
    expect_equal(markovchainFit(as.data.frame(matseq), sanitize = TRUE)$estimate@transitionMatrix, 
                 matrix(c(0, 1, 0, 0, 0, 1, 1/3, 1/3, 1/3), nrow = 3, 
                        byrow = TRUE, dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))

})

### Markovchain Fitting(mle) with sanitize parameter
mle_sequence <- c("a", "b", "b", "a", "a", "a", "b", "b", "b", "a", "a", "b", "a", "a", "b", "c")
mle_fit1 <- markovchainFit(mle_sequence)
mle_fit2 <- markovchainFit(mle_sequence, sanitize = TRUE)

test_that("MarkovchainFit MLE", {
  expect_equal(mle_fit1$estimate@transitionMatrix, 
               matrix(c(0.5, 0.5, 0, 3/7, 3/7, 1/7, 0, 0, 0), nrow = 3, 
                      byrow = TRUE, dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))
  
  expect_equal(mle_fit2$estimate@transitionMatrix, 
               matrix(c(0.5, 0.5, 0, 3/7, 3/7, 1/7, 1/3, 1/3, 1/3), nrow = 3, 
                      byrow = TRUE, dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))
  
  expect_equal(mle_fit1$logLikelihood, mle_fit2$logLikelihood)
  
  expect_equal(mle_fit1$confidenceInterval, mle_fit2$confidenceInterval)
  
  expect_equal(mle_fit2$standardError, mle_fit2$standardError)
})

### Markovchain Fitting(laplace) with sanitize parameter
lap_sequence <- c("a", "b", "b", "a", "a", "a", "b", "b", "b", "a", "a", "b", "a", "a", "b", "c")
lap_fit1 <- markovchainFit(lap_sequence, "laplace")
lap_fit2 <- markovchainFit(lap_sequence, "laplace", sanitize = TRUE)

test_that("Markovchain Laplace", {
  expect_equal(lap_fit1$estimate@transitionMatrix, 
               matrix(c(0.5, 0.5, 0, 3/7, 3/7, 1/7, 0, 0, 0), nrow = 3, 
                      byrow = TRUE, dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))
  
  expect_equal(lap_fit2$estimate@transitionMatrix, 
               matrix(c(0.5, 0.5, 0, 3/7, 3/7, 1/7, 1/3, 1/3, 1/3), nrow = 3, 
                      byrow = TRUE, dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))
  
  expect_equal(lap_fit1$logLikelihood, lap_fit2$logLikelihood)
})

### Markovchain Fitting when some states are not present in the given sequence
mix_seq <- c("a", "b", "b", "a", "a", "a", "b", "b", "b", "a", "a", "b", "a", "a", "b", "c")

mix_fit1 <- markovchainFit(mix_seq, "mle", sanitize = TRUE, possibleStates = c("d"))
mix_fit2 <- markovchainFit(mix_seq, "laplace", sanitize = TRUE, possibleStates = c("d")) 
mix_fit3 <- markovchainFit(mix_seq, "map", sanitize = TRUE, possibleStates = c("d")) 

test_that("Mixture of Markovchain Fitting", {
  expect_equal(mix_fit2$estimate@transitionMatrix, 
               matrix(c(.5, .5, 0, 0, 3/7, 3/7, 1/7, 0,
                        1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4), nrow = 4, byrow = TRUE,
                      dimnames = list(c("a", "b", "c", "d"), c("a", "b", "c", "d"))
                      )
              )
  
  expect_equal(mix_fit1$estimate@transitionMatrix, 
               matrix(c(.5, .5, 0, 0, 3/7, 3/7, 1/7, 0,
                        1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4), nrow = 4, byrow = TRUE,
                      dimnames = list(c("a", "b", "c", "d"), c("a", "b", "c", "d"))
               )
  )
  
  expect_equal(mix_fit3$estimate@transitionMatrix, 
               matrix(c(.5, .5, 0, 0, 3/7, 3/7, 1/7, 0,
                        1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4), nrow = 4, byrow = TRUE,
                      dimnames = list(c("a", "b", "c", "d"), c("a", "b", "c", "d"))
               )
  )
  
})

### Test for createSequenceMatrix
rsequence <- c("a", "b", "b", "a", "a", "a", "b", "b", "b", "a", "a", "b", "a", "a", "b", "c")

test_that("createSequenceMatrix : Permutation of parameters",{
  expect_equal(createSequenceMatrix(rsequence, FALSE, FALSE), 
                   matrix(c(4, 4, 0, 3, 3, 1, 0, 0, 0), nrow = 3, 
                          byrow = TRUE, dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))
  
  expect_equal(createSequenceMatrix(rsequence, FALSE, TRUE), 
               matrix(c(4, 4, 0, 3, 3, 1, 1, 1, 1), nrow = 3, 
                      byrow = TRUE, dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))
  
  expect_equal(createSequenceMatrix(rsequence, TRUE, FALSE), 
               matrix(c(4/8, 4/8, 0, 3/7, 3/7, 1/7, 0, 0, 0), nrow = 3, 
                      byrow = TRUE, dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))
  
  expect_equal(createSequenceMatrix(rsequence, TRUE, TRUE), 
               matrix(c(4/8, 4/8, 0, 3/7, 3/7, 1/7, 1/3, 1/3, 1/3), nrow = 3, 
                      byrow = TRUE, dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))
})


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

pDRes <- c(log(3/2), log(3/2))
names(pDRes) <- c("a", "b")
test_that("priorDistribution must sastisfy", {
  expect_equal(priorDistribution(matrix(c(0.5, 0.5, 0.5, 0.5), 
                                            nrow = 2, 
                                            dimnames = list(c("a", "b"), c("a", "b"))), 
                                     matrix(c(2, 2, 2, 2), 
                                            nrow = 2, 
                                            dimnames = list(c("a", "b"), c("a", "b")))), 
                   pDRes)
})

energyStates <- c("sigma", "sigma_star")
byRow <- TRUE
gen <- matrix(data = c(-3, 3,
                       1, -1), nrow = 2,
              byrow = byRow, dimnames = list(energyStates, energyStates))
molecularCTMC <- new("ctmc", states = energyStates, 
                     byrow = byRow, generator = gen, 
                     name = "Molecular Transition Model")      
test_that("steadyStates must satisfy", {
  expect_identical(steadyStates(molecularCTMC), 
                   matrix(c(1/4, 3/4), nrow = 1, dimnames = list(c(), energyStates)))
})