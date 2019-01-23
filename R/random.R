# Methods to generate random markov chains

# Returns a random stochastic (by row) matrix
randomStochasticMatrix <- function(n) {
  matrix <- matrix(runif(n * n), n)

  result <- apply(matrix, MARGIN = 1, function(row) {
    rowSum <- sum(row)
    
    if (rowSum == 0)
      rep (1/n, n)
    else
      sum / rowSum
  }

  result
}
  
