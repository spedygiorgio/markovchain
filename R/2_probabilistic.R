# given a markovchain object is it possible to reach goal state from 
# a given state

is.accessible <- function(object, from, to) {
  # assume that it is not possible
  out <- FALSE
  
  # names of states
  statesNames <- states(object)
  
  # row number
  fromPos <- which(statesNames == from)
  
  # column number
  toPos<-which(statesNames == to)

  # a logical matrix which will tell the reachability of jth state from ith state
  R <- .commStatesFinderRcpp(object@transitionMatrix)
  
  if(R[fromPos, toPos] == TRUE) {
    out <- TRUE 
  }
  
  return(out)
}

# a markov chain is irreducible if it is composed of only one communicating class
is.irreducible <- function(object) {
  # assuming markovchain chain has more than 1 communicating classes
  out <- FALSE

  # this function will return a list of communicating classes 
  tocheck <- .communicatingClassesRcpp(object)
  
  # only one class implies irreducible markovchain
  if(length(tocheck) == 1) {
    out<-TRUE  
  }
  
  return(out)
}

# what this function will do?
# It calculates the probability to go from given state
# to all other states in k steps
# k varies from 1 to n

firstPassage <- function(object, state, n) {
  P <- object@transitionMatrix
  stateNames <- states(object)
  
  # row number
  i <- which(stateNames == state)

  
  outMatr <- .firstpassageKernelRcpp(P = P, i = i, n = n)
  colnames(outMatr) <- stateNames
  rownames(outMatr) <- 1:n
  return(outMatr)
}

# return a list of communicating classes
communicatingClasses <- function(object) {
  out <- .communicatingClassesRcpp(object)
  return(out)
}

# A communicating class will be a recurrent class if 
# there is no outgoing edge from this class
# Recurrent classes are subset of communicating classes
recurrentClasses <- function(object) {
  out <- .recurrentClassesRcpp(object)
  return(out)
}


#' @title Check if a DTMC is regular
#' 
#' @description Function to check wether a DTCM is regular
#' 
#' @details A regular Markov chain has $A^n$ strictly positive for some n. 
#' So we check: if there is only one eigenvector; if the steadystate vector is striclty positive.
#' 
#' @param object a markovchain object
#' 
#' @return A boolean value
#' 
#' @examples 
#' P=matrix(c(0.5,.25,.25,.5,0,.5,.25,.25,.5),nrow = 3)
#' colnames(P)<-rownames(P)<-c("R","N","S")
#' ciao<-as(P,"markovchain")
#' is.regular(ciao)
#' 
#' @seealso \code{\link{is.irreducible}}


# is.regular<-function(object) {
#   eigenValues<-steadyStates(object = object)
#   minDim<-min(dim(eigenValues))
#   out <- minDim==1 & all(eigenValues>0)
#   return(out)
# }