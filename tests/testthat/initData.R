numInstances <- 100
maxDim <- 100

set.seed(1234567)
randomDims  <- sample(1:100, numInstances, replace = TRUE)

# Get 1:[maxDim] identity by-row-markov-chains
diagonalMCs <- lapply(1:numInstances, function(n) {
  new("markovchain", transitionMatrix = diag(n))
})

# Get [numInstances] random by-row-markov-chains, 
# with dimensions ranging from 1 to 100
MCs <- lapply(randomDims, function(n) {
  randomMarkovChain(n)
})