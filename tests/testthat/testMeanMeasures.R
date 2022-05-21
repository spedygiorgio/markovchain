context("Checking meanFirstPassageTime and meanRecurrenceTime")

# Prepare matrices with known solution

scr <- c("s","c","r")

Pmat <- matrix( c(6,3,1,
                  2,3,5,
                  4,1,5)/10,3,byrow=T)
P <- new("markovchain", states=scr, transitionMatrix=Pmat)

# Analytic solutions
P_r    <- c(s=50,c=30)/11
P_full <- matrix( c( 0,    15/4, 50/11,
                     10/3, 0,    30/11,
                     8/3,  5,    0 ), byrow=T, ncol=3)
rownames(P_full) <- scr
colnames(P_full) <- scr


Poz <- new("markovchain", states=scr,
           transitionMatrix=matrix(c(2,1,1,
                                     2,0,2,
                                     1,1,2)/4, byrow=T, ncol=3))

Poz_full <- matrix( c( 0,  4, 10/3,
                     8/3,  0, 8/3,
                     10/3, 4, 0   ), byrow=T, ncol=3)
rownames(Poz_full) <- scr
colnames(Poz_full) <- scr

test_that("meanFirstPassageTime works for known matrices", {
  expect_equal(meanFirstPassageTime(P,"r"), P_r)
  expect_equal(meanFirstPassageTime(P),     P_full)
  expect_equal(meanFirstPassageTime(Poz),   Poz_full)
})

# Given M = (m_{ij}) where m_{ij} is the mean recurrence time from i to j
# Given P the transition probabilities
# Given C a matrix with all its components as a 1
# Given D a matrix where the diagonal is formed by the recurrence times r_i and
#   the rest of the elements are 0s
#
# It must hold M = PM + C - D (by rows equation)
test_that("meanFirstPassageTime and recurrenceTime hold their characteristic equation", {
  for (mc in allPositiveMCs) {
    P <- mc$transitionMatrix
    M <- meanFirstPassageTime(mc$object)
    C <- markovchain:::ones(ncol(P))
    D <- diag(meanRecurrenceTime(mc$object))

    if (mc$byrow)
      expect_true(all.equal(M, P %*% M + C - D))
    else
      expect_true(all.equal(M, M %*% P + C - D))
  }
})

# Note that meanRecurrenceTimes are the inverse of the steady states elements
# which are not negative,
#
# One steady state:          Other:
#    0                         0
#    0                       u > 0
#    .                       v > 0
#    .                         0
#    .                         .
#   x > 0                      .
#   y > 0                      .
#   z > 0                      .
#    0                         .
#
# So if we invert the mean recurrenceTimes and fill the positions corresponding
# to transient states with 0s, the result should be an eigen vector of the
# transition matrix
#
test_that("We can manufacture an eigen vector with meanRecurrenceTimes", {
  for (mc in allMCs) {
    P <- mc$transitionMatrix
    byrow <- mc$byrow
    times <- mc$meanRecurrenceTime
    states <- mc$states
    inverse <- times ** (-1)
    v <- sapply(states, function(s) {
      ifelse(is.na(inverse[s]), 0, inverse[s])
    })
    v <- unname(v)

    if (byrow)
      result <- as.numeric(v %*% P)
    else
      result <- as.numeric(P %*% v)

    expect_true(all.equal(result, v))
  }
})


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
    #expect_true(all(numVisits >= 0))
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

  M <- zeros(5)
  M[1,1] <- M[5,5] <- 1
  M[2,1] <- M[2,3] <- 1/2
  M[3,2] <- M[3,4] <- 1/2
  M[4,2] <- M[4,5] <- 1/2

  markovChain <- new("markovchain", transitionMatrix = M)

  result <- zeros(5)
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


context("Checking absorptionProbabilities")

test_that("Test mean absorption times for known matrix", {
  # For mcHitting, defined in data-raw/db4Tests.R
  M <- zeros(3)
  result <- 1/5 * matrix(c(4, 1, 3, 2, 2, 3), nrow = 3, byrow = TRUE)
  rownames(result) <- c(2, 3, 4)
  colnames(result) <- c(1, 5)
  
  expect_equal(absorptionProbabilities(mcHitting), result)
  expect_equal(absorptionProbabilities(t(mcHitting)), t(result))
})

# Fs is mean absorption probabilities
# N is the fundamental matrix, (I - Q)^{-1}
# This equation is by rows, need to transpose left part for by columns matrices
test_that("Test that (I - Q) Fs = P[transient, recurrent]", {
  for (mc in allMCs) {
    if (length(mc$transientStates) > 0) {
      recurrent <- mc$recurrentStates
      transient <- mc$transientStates
      byrow <- mc$byrow
      states <- mc$states
      whichRecurrent <- which(states %in% recurrent)
      whichTransient <- which(states %in% transient)
      P <- mc$transitionMatrix
      Fs <- absorptionProbabilities(mc$object)
      Ninv <- diag(length(transient)) - P[whichTransient, whichTransient, drop = FALSE]
      
      if (byrow) {
        expected <- P[whichTransient, whichRecurrent, drop = FALSE]
        expect_equal(Ninv %*% Fs, expected)
      } else {
        expected <- P[whichRecurrent, whichTransient, drop = FALSE]
        expect_equal(Fs %*% Ninv, expected)
      }
    }
  }
})


context("Checking meanAbsorptionTime")


test_that("Mean absorption time for known matrix", {
  result <- c(3, 4, 3)
  names(result) <- c(2, 3, 4)
  
  expect_equal(meanAbsorptionTime(mcDrunkard), result)
  expect_equal(meanAbsorptionTime(t(mcDrunkard)), result)
})


test_that("All mean absorption times are greater or equal than 1", {
  for (mc in allMCs) {
    if (length(mc$transientStates) > 0) {
      expect_true(all(meanAbsorptionTime(mcDrunkard) > 1))
    }
  }
})