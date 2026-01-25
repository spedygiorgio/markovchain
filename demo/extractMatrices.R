#function to extract matrices 
require(markovchain)
extractMatrices <- function(mcObj) {
  
  mcObj <- canonicalForm(object = mcObj)
  #get the indices of transient and absorbing
  transIdx <- which(states(mcObj) %in% transientStates(mcObj))
  absIdx <- which(states(mcObj) %in% absorbingStates(mcObj))
  #get Q, R and I
  Q <- as.matrix(mcObj@transitionMatrix[transIdx,transIdx])
  R <- as.matrix(mcObj@transitionMatrix[transIdx,absIdx])
  I <- as.matrix(mcObj@transitionMatrix[absIdx, absIdx])
  #get fundamental matrix
  N <- solve(eye(size(Q)) - Q)
  #final absorbion probabilities
  NR <- N %*% R
  #return
  out <- list(
    canonicalForm = mcObj,
    Q = Q,
    R = R,
    I = I, 
    N=N,
    NR=NR
  )
  return(out)
}