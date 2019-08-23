context("Checking is.regular")


#https://www.math.ucdavis.edu/~gravner/MAT135B/materials/ch13.pdf
mc1Matrix <- matrix(c(0.5, 0.5, 0,
                      0.5, 0.25, 0.25,
                      0, 1/3, 2/3), 
                    nrow = 3, byrow = TRUE)
mc1 <- as(mc1Matrix, "markovchain")


test_that("Markov chains with strictly positive transition matrices are regular", {
  for (mc in allPositiveMCs) {
    expect_true(mc$regular)
  }
})


test_that("Regularity implies ergodicity", {
  for (mc in allAndPositiveMCs)
    
    if (mc$regular)
      expect_true(mc$irreducible)
})


# Perron–Frobenius theorem: a non negative matrix is primitive 
# (i.e.regular) iff 1 is the maximal unique eigen value
test_that("Regularity iff a single eigen value of modulo |1|", {
  for (mc in allMCs) {
    if (mc$irreducible) {
      # Compute the number of eigen values greater or equal than 1
      eigenValues <- eigen(mc$transitionMatrix, only.values = TRUE)$values
      eigenValues <- sapply(eigenValues, abs)
      maxEigenValues <- sapply(eigenValues, function(e) {isTRUE(all.equal(e, 1)) || e > 1 })
      numMaxEigenValues <- length(which(maxEigenValues))
      
      if (numMaxEigenValues == 1)
        expect_true(mc$regular)
      else
        expect_false(mc$regular)
    }
  }
})


context("Checking canonicForm and is.irreducible")


test_that("Markov chain is irreducible iff there is a single communicating class", {
  
  for (mc in allAndPositiveMCs) {
    states <- mc$states
    commClasses  <- mc$communicatingClasses
    numCommClasses <- length(commClasses)
    irreducible <- mc$irreducible
    
    if (irreducible)
      expect_equal(numCommClasses, 1)
    if (numCommClasses == 1)
      expect_true(irreducible)
  }
})


test_that("If the matrix is irreducible then the canonic form equals the Markov chain", {
  
  for (mc in allAndPositiveMCs) {
    canonic <- mc$canonicForm
    irreducible <- mc$irreducible
    canonicEqual <- canonic == mc$object
    
    if (irreducible)
      expect_true(canonicEqual)
  }
})


test_that("Check known Markov chain is irreducible", {
  expect_true(is.irreducible(mc1))
})
