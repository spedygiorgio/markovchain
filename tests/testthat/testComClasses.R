context("Checking that commClassesKernelRcpp works as expected")


# Not very good in efficiency, but it serves its purpose though
# O(n³) implementation
checkInterchangeability <- function(matrix) {
  # Matrix should be square matrix
  n <- ncol(matrix)
  
  correctCommClasses <- sapply(1:n, function(i) {
    currentRow <- matrix[i, ]
    onesIdx <- which(currentRow)
    
    whichEqual <- sapply(onesIdx, function(j) { matrix[j, ] == currentRow})
    # Is there any row unequal to the row taken as ref
    any(!whichEqual)
  })
  
  any(!correctCommClasses)
}


test_that("Communicating classes matrix is symmetric", {
  
  for (markovChain in MCs) {
    transitionMatrix <- attr(markovChain, "transitionMatrix")
    communicating <- .commClassesKernelRcpp(transitionMatrix)
    C <- communicating$C
    
    expect_equal(C, t(C))
  }
})


test_that("Rows of the same class are interchangeable in a communicating classes matrix", {
  
  for (markovChain in MCs) {
    transitionMatrix <- attr(markovChain, "transitionMatrix")
    communicating <- .commClassesKernelRcpp(transitionMatrix)
    C <- communicating$C
    
    expect_equal(checkInterchangeability(C), TRUE)
  }
})


test_that("Communicating classes of identity matrix of size n are {1, ..., n}", {
  
  for (markovChain in diagonalMCs) {
    transitionMatrix <- attr(markovChain, "transitionMatrix")
    states <- attr(markovChain, "states")
    expected <- as.matrix(apply(transitionMatrix, 1, function(x){ x == 1 }))
    colnames(expected) <- states
    rownames(expected) <- states
    communicating <- .commClassesKernelRcpp(transitionMatrix)
    C <- communicating$C
    
    expect_equal(C, expected)
  }
})


test_that("All clasess are closed for identity matrixes", {
  
  for (markovChain in diagonalMCs) {
    transitionMatrix <- attr(markovChain, "transitionMatrix")
    areClosed <- .commClassesKernelRcpp(transitionMatrix)$v
    
    expect_false(any(!areClosed))
  }
})