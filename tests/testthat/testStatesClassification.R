context("Checking classification of states: recurrentStates, transientStates, absorbingStates")

A <-  structure(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
              0.5, 0.3, 0.3, 0, 0, 0, 0, 0, 0.5, 0.7, 0.7, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0.3, 0.4, 0, 0, 0, 0, 0, 0.4, 0, 0.5, 0, 0, 0, 0, 0,
              0.6, 0.7, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 1), .Dim = c(8L, 8L), .Dimnames = list(
                c("1", "2", "3", "4", "5", "6", "7", "8"), c("1", "2", "3",
                                                             "4", "5", "6", "7", "8")))
mchain <- new("markovchain", transitionMatrix=A)


#https://www.math.ucdavis.edu/~gravner/MAT135B/materials/ch13.pdf
mcMatr1<-matlab::zeros(3)
mcMatr1[1,]<-c(0.5,0.5,0)
mcMatr1[2,]<-c(0.5,0.25,0.25)
mcMatr1[3,]<-c(0,1/3,2/3)
mc1<-as(mcMatr1,"markovchain")


mcMatr2<-matrix(c(0, 0, 1/2, 1/2,1, 0 ,0, 0,0, 1, 0, 0,0, 1, 0, 0),ncol = 4,byrow=TRUE)
mc2<-as(mcMatr2,"markovchain")


mcMatr3<-matrix(c(
0,1,0,0,0,0,
0.4,0.6,0,0,0,0,
0.3,0,0.4,0.2,0.1,0,
0,0,0,0.3,0.7,0,
0,0,0,0.5,0,0.5,
0,0,0,0.3,0,0.7),nrow = 6,byrow=TRUE)

mc3<-as(mcMatr3,"markovchain")
recurrentClasses(mc3)
transientStates(mc3)
mcMatr4<-matlab::zeros(5)
mcMatr4[1:2,1:2]<-0.5*matlab::ones(2)
mcMatr4[5,1]<-1
mcMatr4[3,3]<-1
mcMatr4[4,3:4]<-0.5
mc4<-as(mcMatr4,"markovchain")


test_that("Test recurrent / transient / absorbing states for known Markov chains", {
  expect_equal(recurrentClasses(mc2),list(c("s1","s2","s3","s4")))
  expect_equal(recurrentClasses(mc3),list(c("s1","s2"),c("s4","s5","s6") ))
  expect_equal(transientStates(mc3),"s3")
  expect_equal(recurrentClasses(mc4),list(c("s1","s2"),c("s3")))
  expect_equal(absorbingStates(mc4),"s3")
  expect_equal(transientStates(mc4),c("s4","s5"))
  expect_equal(recurrentClasses(mchain), list(c("3", "4"), c("8")))
  expect_equal(transientStates(mchain), c("1", "2", "5", "6", "7"))
  expect_equal(absorbingStates(mchain), "8")
})


test_that("A state is absorbing iff it is singleton recurrent class", {
  
  for (mc in allMCs) {
    states <- mc$states
    classes <- mc$recurrentClasses
    absorbing <- mc$absorbingStates
    
    expect_true(.testthatAbsorbingAreRecurrentClassRcpp(absorbing, classes))
  }
})


test_that("Recurrent states and transient states are a partition of states", {
  
  for (mc in allMCs) {
    states <- mc$states
    recurrentStates <- mc$recurrentStates
    transientStates <- mc$transientStates
    states <- mc$states
    statesUnion <- sort(unique(append(recurrentStates, transientStates)))
    
    expect_equal(statesUnion, sort(states))
  }
})


test_that("All states are recurrent in a identity Markov chain", {
  
  for (mc in allDiagonalMCs) {
    states <- mc$states
    recurrentStates <- mc$recurrentStates
    expect_true(setequal(recurrentStates, states))
  }
})

test_that("If Markov chain is irreducible then all states are recurrent", {
    
  for (mc in allMCs) {
    states <- mc$states
    recurrent <- mc$recurrentStates
    irreducible <- mc$irreducible
    allRecurrent <- setequal(states, recurrent)
    
    if (irreducible)
      expect_true(allRecurrent)
  }
})


test_that("If there are transient states then Markov chain is not irreducible", {
  
  for (mc in allMCs) {
    states <- mc$states
    transient  <- mc$transientStates
    irreducible <- mc$irreducible
    
    if (length(transient) > 0)
      expect_false(irreducible)
  }
})


context("Checking canonicForm and is.irreducible")


test_that("Markov chain is irreducible iff there is a single communicating class", {
  
  for (mc in allMCs) {
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
  
  for (mc in allMCs) {
    canonic <- mc$canonicForm
    irreducible <- mc$irreducible
    canonicEqual <- canonic == mc$object
    
    if (irreducible)
      expect_true(canonicEqual)
  }
})


context("Checking recurrentClasses method")

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

test_that("Checking recurrent classes for known Markov chains", {
  expect_equal(recurrentClasses(probMc), list(c("a", "c")
                                              , c("b", "g", "i")
                                              , c("f")))
})


test_that("Check know Markov chain is irreducible", {
  expect_true(is.irreducible(mc1))
})


test_that("f(i,j) = 1 for i, j in same recurrent class, hittingProbability(i, k) = 0 for k otherwise", {
  
  for (mc in allMCs) {
    byrow <- mc$byrow
    states <- mc$states
    hitting <- mc$hittingProbabilities
    recurrentClasses <- mc$recurrentClasses
    expect_true(.testthatRecurrentAreHittingRcpp(recurrentClasses, hitting, states, byrow))
  }
})


test_that("Union of recurrentClasses is recurrentStates", {
  
  for (mc in allMCs) {
    recClasses <- mc$recurrentClasses
    recStates <- sort(unlist(recClasses))
    target <- sort(mc$recurrentStates)
    
    expect_equal(recStates, target)
  }
})

test_that("Recurrent classes are disjoint", {
  
  for (mc in allMCs) {
    recClasses <- mc$recurrentClasses
    numRecurrentStates <- sum(sapply(recClasses, function(c){ length(c) }))
    numUnion <- length(unique(unlist(recClasses)))
    
    expect_equal(numRecurrentStates, numUnion)
  }
})

context("Checking transientClasses method")