library(markovchain)

P <- matlab::zeros(10)
P[1, c(1, 3)] <- 1/2;
P[2, 2] <- 1/3; P[2,7] <- 2/3;
P[3, 1] <- 1;
P[4, 5] <- 1;
P[5, c(4, 5, 9)] <- 1/3;
P[6, 6] <- 1;
P[7, 7] <- 1/4; P[7,9] <- 3/4;
P[8, c(3, 4, 8, 10)] <- 1/4;
P[9, 2] <- 1;
P[10, c(2, 5, 10)] <- 1/3;
rownames(P) <- letters[1:10]
colnames(P) <- letters[1:10]
probMc <- new("markovchain", transitionMatrix = P,
             name = "Probability MC")

context("Recurrent Classes")

test_that("Checking recurrent classes for a known matrix", {
  expect_equal(recurrentClasses(probMc), list(c("a", "c")
                                            , c("b", "g", "i")
                                            , c("f")))
})

test_that("f(i,j) = 1 for i, j in same recurrent class, hittingProbability(i, k) = 0 for k otherwise", {
  
  for (markovChain in allMCs) {
    byrow <- markovChain@byrow
    states <- markovChain@states
    hitting <- hittingProbabilities(markovChain)
    recurrentClasses <- recurrentClasses(markovChain)
    expect_true(.testthatRecurrentAreHittingRcpp(recurrentClasses, hitting, states, byrow))
  }
})

test_that("All hitting probabilities are 1 iff the Markov chain is irreducible", {
  
  for (markovChain in allMCs) {
    hitting <- hittingProbabilities(markovChain)
    hittingOne <- .testthatHittingAreOneRcpp(hitting)
    irreducible <- is.irreducible(markovChain)
    
    if (irreducible)
      expect_true(hittingOne)
    if (hittingOne)
      expect_true(irreducible)
  }
})

test_that("Union of recurrentClasses is recurrentStates", {
  
  for (markovChain in allMCs) {
    recClasses <- recurrentClasses(markovChain)
    recStates <- sort(unlist(recClasses))
    target <- sort(recurrentStates(markovChain))

    expect_equal(recStates, target)
  }
})

test_that("Recurrent classes are disjoint", {
  
  for (i in 1:length(allMCs)) {
    markovChain <- allMCs[[i]]
    recClasses <- recurrentClasses(markovChain)
    numRecurrentStates <- sum(sapply(recClasses, function(c){ length(c) }))
    numUnion <- length(unique(unlist(recClasses)))

    expect_equal(numRecurrentStates, numUnion)
  }
})