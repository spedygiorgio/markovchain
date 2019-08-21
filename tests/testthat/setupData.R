numInstances <- 150
# Instances of markov chains with all transitions positive
smallNumInstances <- 20
numByRow <- numInstances / 2
maxDim <- 80
set.seed(1234567)

transposed <- function(markovChains) {
  lapply(markovChains, function(mc) { t(mc) })
}

randomDims  <- sample(1:maxDim, numByRow, replace = TRUE)
# Make positive matrices smaller, since it is costly to make computations on them
randomPositiveDims <- sample(1:20, smallNumInstances, replace = TRUE)

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

.colMCs <- lapply(randomDims, function(s) { randomMarkovChain(s, byrow = FALSE) })

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

subsetAllMCs <- sample(allMCs, smallNumInstances, replace = FALSE)

steadyStatesMCs <- lapply(knownSteadyStatesMCs, markovchain:::precomputeData)

allDiagonalMCs <- lapply(.allDiagonalMCs, markovchain:::precomputeData)

allAndDiagonalMCs <- append(allMCs, allDiagonalMCs)

allPositiveMCs <- lapply(.allPositiveMCs, function(mc) {
  list(object = mc,
       byrow = mc@byrow,
       states = mc@states,
       communicatingClasses = communicatingClasses(mc),
       canonicForm = canonicForm(mc),
       transitionMatrix = mc@transitionMatrix,
       regular = is.regular(mc),
       irreducible = is.irreducible(mc)
  )
})

allAndPositiveMCs <- append(allMCs, allPositiveMCs)