#library(markovchain)

test_that("Hitting probabilities of identity markov chain is identity", {
  
  for (markovChain in diagonalMCs) {
    transitionMatrix <- attr(markovChain, "transitionMatrix")
    numStates <- nrow(transitionMatrix)
    expect_equal(hittingProbabilities(markovChain), diag(numStates))
  }
})

test_that("Hitting probabilities of known markov chain", {
  
  M <- matlab::zeros(5, 5)
  M[1,1] <- M[5,5] <- 1
  M[2,1] <- M[2,3] <- 1/2
  M[3,2] <- M[3,4] <- 1/2
  M[4,2] <- M[4,5] <- 1/2
  
  markovChain <- new("markovchain", transitionMatrix = M)
  hittingProbabilities(markovChain)
  
  result <- matlab::zeros(5, 5)
  result[1,1] <- result[5,5] <- 1
  result[2,1] <- 4/5
  result[3,1] <- 3/5
  result[4,1] <- 2/5
  result[2,2] <- 3/8
  result[3,2] <- 3/4
  result[4,2] <- 1/2
  result[2,3] <- 1/2
  result[3,3] <- 3/8
  result[4,3] <- 1/4
  result[2,4] <- 1/3
  result[3,4] <- 2/3
  result[4,4] <- 1/6
  result[2,5] <- 1/5
  result[3,5] <- 2/5
  result[4,5] <- 3/5
  
  expect_equal(hittingProbabilities(markovChain), result)
})