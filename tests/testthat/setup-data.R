numInstances <- 150
numByRow <- numInstances / 2
maxDim <- 100
set.seed(1234567)

transposed <- function(markovChains) {
  lapply(markovChains, function(mc) { t(mc) })
}

randomDims  <- sample(1:maxDim, numByRow, replace = TRUE)

# Get 1:[maxDim] identity by-row-markov-chains
diagonalMCs <- lapply(1:numByRow, function(n) {
  new("markovchain", transitionMatrix = diag(n))
})

colDiagonalMCs <- transposed(diagonalMCs)

# Append by-columns MarkovChains
allDiagonalMCs <- append(diagonalMCs, colDiagonalMCs)

# Get [numByRow] random by-row-markov-chains, 
# with dimensions ranging from 1 to 100
MCs <- lapply(randomDims, function(n) {
  randomMarkovChain(n)
})

colMCs <- transposed(MCs)

# Append by-columns MarkovChains
allMCs <- append(MCs, colMCs)