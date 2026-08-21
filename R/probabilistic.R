# given a markovchain object is it possible to reach goal state from 
# a given state

#' @name is.accessible
#' @title Verify if a state j is reachable from state i.
#' @description This function verifies if a state is reachable from another, i.e., 
#'              if there exists a path that leads to state j leaving from state i with 
#'              positive probability
#'              
#' @param object A \code{markovchain} object.
#' @param from The name of state "i" (beginning state).
#' @param to The name of state "j" (ending state).
#' 
#' @details It wraps an internal function named \code{reachabilityMatrix}.
#' @return A boolean value.
#' 
#' @references James Montgomery, University of Madison
#' 
#' @author Giorgio Spedicato, Ignacio Cordón
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
#' 
#' @exportMethod is.accessible
setGeneric("is.accessible", function(object, from, to) standardGeneric("is.accessible"))

setMethod("is.accessible", c("markovchain", "character", "character"), 
  function(object, from, to) {
    # O(n²) procedure to see if to state is reachable starting at from state
    return(.isAccessibleRcpp(object, from, to))
  }
)

setMethod("is.accessible", c("markovchain", "missing", "missing"), 
  function(object, from, to) {
    .reachabilityMatrixRcpp(object)
  }
)

# a markov chain is irreducible if it is composed of only one communicating class

#' @name is.irreducible
#' @title Function to check if a Markov chain is irreducible (i.e. ergodic)
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
#' 
#' @exportMethod is.irreducible
setGeneric("is.irreducible", function(object) standardGeneric("is.irreducible"))

setMethod("is.irreducible", "markovchain", function(object) {
  .isIrreducibleRcpp(object)
})

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

#' function to calculate first passage probabilities
#' 
#' @description The function calculates first passage probability for a subset of
#' states given an initial state.
#' 
#' @param object a markovchain-class object
#' @param state intital state of the process (charactervector)
#' @param set set of states A, first passage of which is to be calculated
#' @param n Number of rows on which compute the distribution
#' 
#' @return A vector of size n showing the first time proabilities
#' @references
#' Renaldo Feres, Notes for Math 450 Matlab listings for Markov chains;
#' MIT OCW, course - 6.262, Discrete Stochastic Processes, course-notes, chap -05
#' 
#' @author Vandit Jain
#' 
#' @seealso \code{\link{firstPassage}}
#' @examples 
#' statesNames <- c("a", "b", "c")
#' markovB <- new("markovchain", states = statesNames, transitionMatrix =
#' matrix(c(0.2, 0.5, 0.3,
#'          0, 1, 0,
#'          0.1, 0.8, 0.1), nrow = 3, byrow = TRUE,
#'        dimnames = list(statesNames, statesNames)
#' ))
#' firstPassageMultiple(markovB,"a",c("b","c"),4)  
#' 
#' @export 
firstPassageMultiple <- function(object,state,set, n){
  
  # gets the transition matrix
  P <- object@transitionMatrix
  
  # character vector of states of the markovchain
  stateNames <- states(object)
  
  k <- -1
  k <- which(stateNames == state)
  if(k==-1)
    stop("please provide a valid initial state")
  
  # gets the set in numeric vector
  setno <- rep(0,length(set))
  for(i in 1:length(set))
  {
    setno[i] = which(set[i] == stateNames)
    if(setno[i] == 0)
      stop("please provide proper set of states")
  }
  
  # calls Rcpp implementation
  outMatr <- .firstPassageMultipleRCpp(P,k,setno,n)
  
  #sets column and row names of output
  colnames(outMatr) <- "set"
  rownames(outMatr) <- 1:n
  return(outMatr)
}

#' @name communicatingClasses
#' @rdname structuralAnalysis
#' @aliases transientStates
#' @aliases recurrentStates
#' @aliases absorbingStates
#' @aliases communicatingClasses
#' @aliases transientClasses
#' @aliases recurrentClasses
#' @title Various function to perform structural analysis of DTMC
#' @description These functions return absorbing and transient states of the \code{markovchain} objects.
#' 
#' @param object A \code{markovchain} object.
#' 
#' @return
#' \describe{
#'   \item{\code{period}}{returns a integer number corresponding to the periodicity of the Markov 
#'     chain (if it is irreducible)}
#'   \item{\code{absorbingStates}}{returns a character vector with the names of the absorbing 
#'     states in the Markov chain}
#'   \item{\code{communicatingClasses}}{returns a list in which each slot contains the names of
#'     the states that are in that communicating class}
#'   \item{\code{recurrentClasses}}{analogously to \code{communicatingClasses}, but with 
#'     recurrent classes}
#'   \item{\code{transientClasses}}{analogously to \code{communicatingClasses}, but with 
#'     transient classes}
#'   \item{\code{transientStates}}{returns a character vector with all the transient states
#'     for the Markov chain}
#'   \item{\code{recurrentStates}}{returns a character vector with all the recurrent states 
#'     for the Markov chain}
#'   \item{\code{canonicForm}}{returns the Markov chain reordered by a permutation of states 
#'     so that we have blocks submatrices for each of the recurrent classes and a collection 
#'     of rows in the end for the transient states}
#' }
#' 
#' @references Feres, Matlab listing for markov chain.
#' 
#' @author Giorgio Alfredo Spedicato, Ignacio Cordón
#' 
#' @seealso \code{\linkS4class{markovchain}}
#' 
#' @examples 
#' statesNames <- c("a", "b", "c")
#' mc <- new("markovchain", states = statesNames, transitionMatrix =
#'           matrix(c(0.2, 0.5, 0.3,
#'                    0,   1,   0,
#'                    0.1, 0.8, 0.1), nrow = 3, byrow = TRUE,
#'                  dimnames = list(statesNames, statesNames))
#'          )
#' 
#' communicatingClasses(mc)
#' recurrentClasses(mc)
#' recurrentClasses(mc)
#' absorbingStates(mc)
#' transientStates(mc)
#' recurrentStates(mc)
#' canonicForm(mc)
#' 
#' # periodicity analysis
#' A <- matrix(c(0, 1, 0, 0, 0.5, 0, 0.5, 0, 0, 0.5, 0, 0.5, 0, 0, 1, 0), 
#'             nrow = 4, ncol = 4, byrow = TRUE)
#' mcA <- new("markovchain", states = c("a", "b", "c", "d"), 
#'           transitionMatrix = A,
#'           name = "A")
#'
#' is.irreducible(mcA) #true
#' period(mcA) #2
#'
#' # periodicity analysis
#' B <- matrix(c(0, 0, 1/2, 1/4, 1/4, 0, 0,
#'                    0, 0, 1/3, 0, 2/3, 0, 0,
#'                    0, 0, 0, 0, 0, 1/3, 2/3,
#'                    0, 0, 0, 0, 0, 1/2, 1/2,
#'                    0, 0, 0, 0, 0, 3/4, 1/4,
#'                    1/2, 1/2, 0, 0, 0, 0, 0,
#'                    1/4, 3/4, 0, 0, 0, 0, 0), byrow = TRUE, ncol = 7)
#' mcB <- new("markovchain", transitionMatrix = B)
#' period(mcB)
#' 
#' @exportMethod communicatingClasses
setGeneric("communicatingClasses", function(object) standardGeneric("communicatingClasses"))

setMethod("communicatingClasses", "markovchain", function(object) {
  return(.communicatingClassesRcpp(object))
})

# A communicating class will be a recurrent class if 
# there is no outgoing edge from this class
# Recurrent classes are subset of communicating classes
#' @rdname structuralAnalysis
#' 
#' @exportMethod recurrentClasses
setGeneric("recurrentClasses", function(object) standardGeneric("recurrentClasses"))

setMethod("recurrentClasses", "markovchain", function(object) {
  return(.recurrentClassesRcpp(object))
})

# A communicating class will be a transient class iff
# there is an outgoing edge from this class to an state
# outside of the class
# Transient classes are subset of communicating classes
#' @rdname structuralAnalysis
#' 
#' @exportMethod transientClasses
setGeneric("transientClasses", function(object) standardGeneric("transientClasses"))

setMethod("transientClasses", "markovchain", function(object) {
  return(.transientClassesRcpp(object))
})
