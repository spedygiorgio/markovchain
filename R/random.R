# Methods to generate random markov chains

normalizeMatrix <- function(matrix, byrow = TRUE) {
  margin <- ifelse (byrow, 1, 2)
  n <- nrow(matrix)
  
  result <- apply(matrix, MARGIN = margin, function(row) {
    rowSum <- sum(row)
    
    if (rowSum == 0)
      rep (1/n, n)
    else
      row / rowSum
  })
  
  # If we want the result by rows, we have to transpose the matrix,
  # since the apply method with margin = 1 (over rows) returns the result
  # by columns
  if (byrow)
    t(result)
  else
    result
}


# Returns a random stochastic matrix
randomStochasticMatrix <- function(n, byrow = TRUE) {
  numRandom <- n * n
  randomNums <- runif(numRandom)
  zeroProbs  <- runif(n)
  
  # Adjust to have 0s in each row
  result <- sapply(1:n, function(i) {
    zeroProb <- zeroProbs[i]
    remainProb <- (1 - zeroProb) / numRandom
    probs <- c(zeroProb, rep(remainProb, numRandom))
    rowEntries <- sample(c(0, randomNums), n, prob = probs, replace = TRUE)
    
    rowEntries
  })
  
  if (byrow)
    result <- t(result)
  
  result <- normalizeMatrix(result, byrow)
  
  result
}


randomMarkovChain <- function(n) {
  matrix <- randomStochasticMatrix(n)
  new("markovchain", transitionMatrix = matrix)
}