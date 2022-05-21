context("Checking hittingProbabilities method")

test_that("Hitting probabilities of identity markov chain is identity", {
  
  for (mc in allDiagonalMCs) {
    states <- mc$states
    numStates <- length(states)
    result <- diag(numStates)
    rownames(result) <- states
    colnames(result) <- states
    hittingProbabilities <- mc$hittingProbabilities
    
    expect_equal(hittingProbabilities, result)
  }
})


test_that("Hitting probabilities hold their characteristic system and are non negative", {
  # Check that the following recurrence holds,
  # naming p = probs, f = hitting, it checks:
  #
  # f(i, j) = p(i, j) + ∑_{k ≠ j} p(i, k) f(k, j)
  
  for (mc in allMCs) {
    probs <- mc$transitionMatrix
    byrow <- mc$byrow
    hitting <- mc$hittingProbabilities
    #expect_true(all(hitting >= 0))
    expect_true(.testthatAreHittingRcpp(probs, hitting, byrow))
  }
})


test_that("All hitting probabilities are 1 iff the Markov chain is irreducible", {
  
  for (mc in allMCs) {
    hitting <- mc$hittingProbabilities
    hittingOne <- .testthatHittingAreOneRcpp(hitting)
    irreducible <- mc$irreducible
    
    if (irreducible)
      expect_true(hittingOne)
    if (hittingOne)
      expect_true(irreducible)
  }
})


# Test with a matrix with known hitting probabilities
# Taken from the book Procesos Estocásticos, Ricardo Vélez & Tomás Prieto
test_that("Tests hitting probabilities for a known markov chain", {
  # For mcHitting defined in data-raw/db4Tests.R
  
  result <- markovchain:::zeros(5, 5)
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
  rownames(result) <- mcHitting@states
  colnames(result) <- mcHitting@states
  
  expect_equal(hittingProbabilities(mcHitting), result)
  expect_equal(hittingProbabilities(t(mcHitting)), t(result))
})