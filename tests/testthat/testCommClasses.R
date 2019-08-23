context("Checking .commClassesKernelRcpp method")


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
  
  for (mc in allMCs) {
    transitionMatrix <- mc$transitionMatrix
    communicating <- .commClassesKernelRcpp(transitionMatrix)
    C <- communicating$classes
    
    expect_equal(C, t(C))
  }
})


test_that("Rows of the same class are interchangeable in a communicating classes matrix", {
  
  for (mc in allMCs) {
    transitionMatrix <- mc$transitionMatrix
    communicating <- .commClassesKernelRcpp(transitionMatrix)
    C <- communicating$classes
    
    expect_equal(checkInterchangeability(C), TRUE)
  }
})


test_that("Communicating classes of identity matrix of size n are {1, ..., n}", {
  
  for (mc in allDiagonalMCs) {
    transitionMatrix <- mc$transitionMatrix
    states <- mc$states
    expected <- as.matrix(apply(transitionMatrix, 1, function(x){ x == 1 }))
    colnames(expected) <- states
    rownames(expected) <- states
    communicating <- .commClassesKernelRcpp(transitionMatrix)
    C <- communicating$classes
    
    expect_equal(C, expected)
  }
})


test_that("All clasess are closed for identity matrixes", {
  
  for (mc in allDiagonalMCs) {
    transitionMatrix <- mc$transitionMatrix
    areClosed <- .commClassesKernelRcpp(transitionMatrix)$closed
    
    expect_true(all(areClosed))
  }
})


test_that("Communicating class matrix is correct", {
  
  for (mc in allMCs) {
    # P
    transitionMatrix <- mc$transitionMatrix
    n <- ncol(transitionMatrix)
    # The communicating matrix has a 1 in an entry (i,j) iff
    # P'^{n - 1} has a positive number in its entries (i,j) and (j,i)
    # When we say P' we refer to making i always communicate with itself
    p_n <- (transitionMatrix + diag(n)) %^% (n - 1) > 0
    commClasses <- .commClassesKernelRcpp(transitionMatrix)$classes
    # Correct the diagonal to be always positive 
    # (i always communicates with itself)
    expectedCommMatrix <- (p_n * t(p_n)) > 0
    
    expect_true(all(commClasses == expectedCommMatrix))
  }
})


context("Checking communicatingClasses method")


test_that("Communicating classes are a partition of states", {
  
  for (mc in allAndDiagonalMCs) {
    states <- mc$states
    commClasses <- mc$communicatingClasses
    expect_true(.testthatIsPartitionRcpp(commClasses, states))
  }
})
