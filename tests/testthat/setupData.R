numInstances <- 150
# Instances of markov chains with all transitions positive
numPositiveInstances <- 20
numByRow <- numInstances / 2
maxDim <- 100
set.seed(1234567)

transposed <- function(markovChains) {
  lapply(markovChains, function(mc) { t(mc) })
}

randomDims  <- sample(1:maxDim, numByRow, replace = TRUE)
randomPositiveDims <- sample(1:maxDim, numPositiveInstances, replace = TRUE)

# Get 1:[maxDim] identity by-row-markov-chains
.diagonalMCs <- lapply(1:numByRow, function(n) {
  new("markovchain", transitionMatrix = diag(n))
})

.colDiagonalMCs <- transposed(.diagonalMCs)

# Append by-columns MarkovChains
.allDiagonalMCs <- append(.diagonalMCs, .colDiagonalMCs)

# Get [numByRow] random by-row-markov-chains, 
# with dimensions ranging from 1 to 100
.MCs <- lapply(randomDims, randomMarkovChain)

.colMCs <- transposed(.MCs)

# Append by-columns MarkovChains
.allMCs <- append(.MCs, .colMCs)

mcsIndexes <- seq_along(.allMCs)
diagonalIndexes <- seq_along(.allDiagonalMCs)

# Markov chains with transition matrix P > 0
.positiveMCs <- lapply(randomPositiveDims, function(n) { randomMarkovChain(n, zeroProb = 0) })
.colPositiveMCs <- transposed(.positiveMCs)
.allPositiveMCs <- append(.positiveMCs, .colPositiveMCs)

#################################################################
# Classes and states pre-computed data
#################################################################

allMCs <- lapply(.allMCs, markovchain:::precomputeData)
steadyStatesMCs <- lapply(knownSteadyStatesMCs, markovchain:::precomputeData)
allDiagonalMCs <- lapply(.allDiagonalMCs, markovchain:::precomputeData)
allAndDiagonalMCs <- append(allMCs, allDiagonalMCs)
allPositiveMCs <- lapply(.positiveMCs, function(mc) {
  list(object = mc,
       regular = is.regular(mc))
})