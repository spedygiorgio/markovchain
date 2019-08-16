context("Checking meanNumVisits method")

test_that("Mean number visits of identity markov chain is identity * Inf", {
  
  for (mc in allDiagonalMCs) {
    states <- mc$states
    numStates <- length(states)
    result <- diag(numStates)
    diag(result) <- Inf
    rownames(result) <- states
    colnames(result) <- states
    meanNumVisits <- mc$meanNumVisits
    
    expect_equal(meanNumVisits, result)
  }
})


test_that("Mean number of visits hold their characteristic system and are non negative", {
  # Check that the following recurrence holds,
  # naming p = probs, f = hitting, E = mean number of visits, it checks:
  #
  # E(i, j) = p(i, j) / (1 - f(j, j)) + ∑_{k ≠ j} p(i, k) E(k, j)
  
  for (mc in allMCs) {
    probs <- mc$transitionMatrix
    byrow <- mc$byrow
    hitting <- mc$hittingProbabilities
    numVisits <- mc$meanNumVisits
    expect_true(all(numVisits >= 0))
    expect_true(.testthatAreMeanNumVisitsRcpp(probs, numVisits, hitting, byrow))
  }
})



test_that("All mean number of visits are ∞ iff the Markov chain is irreducible", {
  
  for (mc in allMCs) {
    meanNumVisits <- mc$meanNumVisits
    numVisitsInf <- all(meanNumVisits == Inf)
    irreducible <- mc$irreducible
    
    if (irreducible)
      expect_true(numVisitsInf)
    if (numVisitsInf)
      expect_true(irreducible)
  }
})


# Test mean number of visits with a known matrix
# Taken from the book Procesos Estocásticos, Ricardo Vélez & Tomás Prieto
test_that("Tests mean number of visits for a known markov chain", {
  
  M <- matlab::zeros(5, 5)
  M[1,1] <- M[5,5] <- 1
  M[2,1] <- M[2,3] <- 1/2
  M[3,2] <- M[3,4] <- 1/2
  M[4,2] <- M[4,5] <- 1/2
  
  markovChain <- new("markovchain", transitionMatrix = M)

  result <- matlab::zeros(5, 5)
  result[1:4, 1] <- Inf
  result[2:5, 5] <- Inf
  result[1, 2:5] <- 0
  result[5, 1:4] <- 0
  result[2,2] <- result[3,3] <- 3/5
  result[2,3] <- result[4,2] <- result[3,4] <- 4/5
  result[2,4] <- result[4,3] <- 2/5
  result[3,2] <- 6/5
  result[4,4] <- 1/5
  rownames(result) <- markovChain@states
  colnames(result) <- markovChain@states
  
  expect_equal(meanNumVisits(markovChain), result)
  expect_equal(meanNumVisits(t(markovChain)), t(result))
})