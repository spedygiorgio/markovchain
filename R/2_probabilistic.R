# given a markovchain object is it possible to reach goal state from 
# a given state

#' @name is.accessible
#' @title Verify if a state j is reachable from state i.
#' @description This function verifies if a state is reachable from another, i.e., 
#'              if exists a path that leads to state j leaving from state i with 
#'              positive probability
#'              
#' @param object A \code{markovchain} object.
#' @param from The name of state "i" (beginning state).
#' @param to The name of state "j" (ending state).
#' 
#' @details It wraps an internal function named \code{.commStatesFinder}.
#' @return A boolean value.
#' 
#' @references James Montgomery, University of Madison
#' 
#' @author Giorgio Spedicato
#' @seealso \code{is.irreducible}
#' 
#' @examples 
#' statesNames <- c("a", "b", "c")
#' markovB <- new("markovchain", states = statesNames, 
#'                transitionMatrix = matrix(c(0.2, 0.5, 0.3,
#'                                              0,   1,   0,
#'                                            0.1, 0.8, 0.1), nrow = 3, byrow = TRUE, 
#'                                          dimnames = list(statesNames, statesNames)
#'                                         )
#'                )
#' is.accessible(markovB, "a", "c")
#' @export

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

#' @name is.irreducible
#' @title Function to check if a Markov chain is irreducible
#' @description This function verifies whether a \code{markovchain} object transition matrix 
#'              is composed by only one communicating class.
#' @param object A \code{markovchain} object
#' 
#' @details It is based on \code{.communicatingClasses} internal function.
#' @return A boolean values.
#' 
#' @references Feres, Matlab listings for Markov Chains.
#' @author Giorgio Spedicato
#' 
#' @seealso \code{\link{summary}}
#' 
#' @examples 
#' statesNames <- c("a", "b")
#' mcA <- new("markovchain", transitionMatrix = matrix(c(0.7,0.3,0.1,0.9),
#'                                              byrow = TRUE, nrow = 2, 
#'                                              dimnames = list(statesNames, statesNames)
#'            ))
#' is.irreducible(mcA)
#' @export

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

#' @name firstPassage
#' @title First passage across states
#' @description This function compute the first passage probability in states
#' 
#' @param object A \code{markovchain} object
#' @param state Initial state
#' @param n Number of rows on which compute the distribution
#' 
#' @details Based on Feres' Matlab listings
#' @return A matrix of size 1:n x number of states showing the probability of the 
#'         first time of passage in states to be exactly the number in the row.
#'
#' @references Renaldo Feres, Notes for Math 450 Matlab listings for Markov chains
#' 
#' @author Giorgio Spedicato
#' @seealso \code{\link{conditionalDistribution}}
#' 
#' @examples 
#' simpleMc <- new("markovchain", states = c("a", "b"),
#'                  transitionMatrix = matrix(c(0.4, 0.6, .3, .7), 
#'                                     nrow = 2, byrow = TRUE))
#' firstPassage(simpleMc, "b", 20)
#'
#' @export

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

#' @name communicatingClasses
#' @title Various function to perform structural analysis of DTMC
#' @description These functions return absorbing and transient states of the \code{markovchain} objects.
#' 
#' @param object A \code{markovchain} object.
#' 
#' @return vector, matrix or list
#' 
#' @references Feres, Matlab listing for markov chain.
#' 
#' @author Giorgio Alfredo Spedicato
#' 
#' @seealso \code{\linkS4class{markovchain}}
#' 
#' @examples 
#' statesNames <- c("a", "b", "c")
#' markovB <- new("markovchain", states = statesNames, transitionMatrix =
#'                    matrix(c(0.2, 0.5, 0.3,
#'                               0,   1,   0,
#'                             0.1, 0.8, 0.1), nrow = 3, byrow = TRUE, 
#'                             dimnames = list(statesNames, statesNames)
#'               ))
#'               
#' communicatingClasses(markovB)               
#' recurrentClasses(markovB)
#' absorbingStates(markovB)
#' transientStates(markovB)
#' canonicForm(markovB)
#' 
#' # periodicity analysis : 1
#' E <- matrix(c(0, 1, 0, 0, 0.5, 0, 0.5, 0, 0, 0.5, 0, 0.5, 0, 0, 1, 0), 
#'             nrow = 4, ncol = 4, byrow = TRUE)
#' mcE <- new("markovchain", states = c("a", "b", "c", "d"), 
#'           transitionMatrix = E, 
#'           name = "E")
#'
#' is.irreducible(mcE) #true
#' period(mcE) #2
#'
#' # periodicity analysis : 2
#' myMatr <- matrix(c(0, 0, 1/2, 1/4, 1/4, 0, 0,
#'                    0, 0, 1/3, 0, 2/3, 0, 0,
#'                    0, 0, 0, 0, 0, 1/3, 2/3,
#'                    0, 0, 0, 0, 0, 1/2, 1/2,
#'                    0, 0, 0, 0, 0, 3/4, 1/4,
#'                    1/2, 1/2, 0, 0, 0, 0, 0,
#'                    1/4, 3/4, 0, 0, 0, 0, 0), byrow = TRUE, ncol = 7)
#' myMc <- new("markovchain", transitionMatrix = myMatr)
#' period(myMc)
#' 
#' @rdname absorbingStates
#' @export

communicatingClasses <- function(object) {
  out <- .communicatingClassesRcpp(object)
  return(out)
}

# A communicating class will be a recurrent class if 
# there is no outgoing edge from this class
# Recurrent classes are subset of communicating classes

#' @rdname absorbingStates
#' 
#' @export

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