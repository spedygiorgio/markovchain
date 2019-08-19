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
    C <- matlab::ones(ncol(P))
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
