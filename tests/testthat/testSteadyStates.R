context("Checking steadyStates method")


steadyStatesMCs <- append(allAndDiagonalMCs, steadyStatesMCs)


test_that("Num of steady states is the same as num of recurrent classes", {
  
  for (mc in steadyStatesMCs) {
    byrow <- mc$byrow
    steady <- mc$steadyStates
    numSteadyStates <- ifelse(byrow, nrow(steady), ncol(steady))
    numRecurrentClasses <- length(mc$recurrentClasses)
    
    expect_equal(numSteadyStates, numRecurrentClasses)
  }
})


test_that("Steady states are prob vectors", {
  
  for (mc in steadyStatesMCs) {
    byrow <- mc$byrow
    steady <- mc$steadyStates
    margin <- ifelse(byrow, 1, 2)
    steadyAreProbVectors <- all(apply(steady, MARGIN = margin, .isProbabilityVector))
    
    expect_true(steadyAreProbVectors)
  }
})


test_that("Steady states are linearly independent vectors", {
  
  for (mc in steadyStatesMCs) {
    byrow <- mc$byrow
    steady <- mc$steadyStates
    rank <- rankMatrix(steady)[[1]]
    
    expect_equal(rank, min(nrow(steady), ncol(steady)))
  }
})


test_that("Steady states v are eigen vectors, i.e. vP = v (by rows) or Pv = v (by cols)", {
  
  for (mc in steadyStatesMCs) {
    byrow <- mc$byrow
    steady <- mc$steadyStates
    P <- mc$transitionMatrix
    margin <- ifelse(byrow, 1, 2)
    areEigenVectors <- apply(steady, MARGIN = margin, function(v) {
      v <- as.numeric(v)
      
      if (byrow)
        result <- as.numeric(v %*% P)
      else
        result <- as.numeric(P %*% v)
      
      all.equal(result, v)
    })
    
    expect_true(all(areEigenVectors))
  }
})

