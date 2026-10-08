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
  if (!is(object, "markovchain")) {
    stop("object must be a markovchain object")
  }
  stateNames <- states(object)
  if (length(state) != 1L || is.na(state) || !state %in% stateNames) {
    stop("state must identify exactly one state of the Markov chain")
  }
  if (length(n) != 1L || is.na(n) || !is.finite(n) ||
      n < 1 || n != floor(n)) {
    stop("n must be a positive integer")
  }

  outMatr <- .firstpassageKernelRcpp(
    P = .rowStochasticMatrix(object),
    i = match(state, stateNames),
    n = as.integer(n)
  )
  colnames(outMatr) <- stateNames
  rownames(outMatr) <- seq_len(nrow(outMatr))
  outMatr
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
firstPassageMultiple <- function(object, state, set, n) {
  if (!is(object, "markovchain")) {
    stop("object must be a markovchain object")
  }
  stateNames <- states(object)
  if (length(state) != 1L || is.na(state) || !state %in% stateNames) {
    stop("state must identify exactly one initial state")
  }
  if (!is.character(set) || length(set) < 1L || anyNA(set) ||
      any(!set %in% stateNames)) {
    stop("set must contain valid state names")
  }
  if (length(n) != 1L || is.na(n) || !is.finite(n) ||
      n < 1 || n != floor(n)) {
    stop("n must be a positive integer")
  }

  out <- .firstPassageMultipleRCpp(
    .rowStochasticMatrix(object),
    match(state, stateNames),
    match(unique(set), stateNames),
    as.integer(n)
  )
  matrix(out, ncol = 1L,
         dimnames = list(seq_along(out), "set"))
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


#' @rdname structuralAnalysis
#' 
#' @exportMethod transientStates
setGeneric("transientStates", function(object) standardGeneric("transientStates"))


setMethod("transientStates", "markovchain", function(object) {
    .transientStatesRcpp(object)
  }
)


#' @rdname structuralAnalysis
#' 
#' @exportMethod recurrentStates
setGeneric("recurrentStates", function(object) standardGeneric("recurrentStates"))


setMethod("recurrentStates", "markovchain", function(object) {
    .recurrentStatesRcpp(object)
  }
)

# generic function to extract absorbing states

#' @rdname structuralAnalysis
#' 
#' @exportMethod absorbingStates
setGeneric("absorbingStates", function(object) standardGeneric("absorbingStates"))

setMethod("absorbingStates", "markovchain", function(object) {
    .absorbingStatesRcpp(object)
  }
)


#' @rdname structuralAnalysis
#' 
#' @exportMethod canonicForm
setGeneric("canonicForm", function(object) standardGeneric("canonicForm"))

setMethod("canonicForm", "markovchain", function(object) {
    .canonicFormRcpp(object)
  }
)


#' @title Calculates committor of a markovchain object with respect to set A, B
#' 
#' @description Returns the probability of hitting states rom set A before set B 
#' with different initial states
#' 
#' @usage committorAB(object,A,B,p)
#' 
#' @param object a markovchain class object
#' @param A a set of states
#' @param B a set of states
#' @param p initial state (default value : 1)
#' 
#' @details The function solves a system of linear equations to calculate probaility that the process hits
#' a state from set A before any state from set B
#' 
#' @return Return a vector of probabilities in case initial state is not provided else returns a number
#' 
#' @examples 
#' transMatr <- matrix(c(0,0,0,1,0.5,
#'                       0.5,0,0,0,0,
#'                       0.5,0,0,0,0,
#'                       0,0.2,0.4,0,0,
#'                       0,0.8,0.6,0,0.5),
#'                       nrow = 5)
#' object <- new("markovchain", states=c("a","b","c","d","e"),transitionMatrix=transMatr)
#' committorAB(object,c(5),c(3))
#' 
#' @export
committorAB <- function(object, A, B, p = 1) {
  if (!is(object, "markovchain")) {
    stop("object must be a markovchain object")
  }

  nstates <- length(object@states)
  valid_indices <- function(x) {
    is.numeric(x) && length(x) > 0L && !anyNA(x) &&
      all(is.finite(x)) && all(x == floor(x)) &&
      all(x >= 1L & x <= nstates)
  }
  if (!valid_indices(A)) stop("please provide a valid set A")
  if (!valid_indices(B)) stop("please provide a valid set B")
  A <- unique(as.integer(A))
  B <- unique(as.integer(B))
  if (length(intersect(A, B)) > 0L) {
    stop("sets A and B must be disjoint")
  }

  return_all <- missing(p)
  if (!return_all &&
      (length(p) != 1L || is.na(p) || !is.finite(p) ||
       p != floor(p) || p < 1L || p > nstates)) {
    stop("please provide a valid initial state")
  }

  coefficient <- .rowStochasticMatrix(object) - diag(nstates)
  coefficient[A, ] <- 0
  coefficient[cbind(A, A)] <- 1
  coefficient[B, ] <- 0
  coefficient[cbind(B, B)] <- 1
  rhs <- numeric(nstates)
  rhs[A] <- 1
  out <- solve(coefficient, rhs)

  if (return_all) out else out[as.integer(p)]
}


#' Expected Rewards for a markovchain
#' 
#' @description Given a markovchain object and reward values for every state,
#' function calculates expected reward value after n steps.
#' 
#' @usage expectedRewards(markovchain,n,rewards)
#' 
#' @param markovchain the markovchain-class object
#' @param n no of steps of the process
#' @param rewards vector depicting rewards coressponding to states
#' 
#' @details the function uses a dynamic programming approach to solve a 
#' recursive equation described in reference.
#' 
#' @return
#' returns a vector of expected rewards for different initial states
#' 
#' @author Vandit Jain
#' 
#' @references Stochastic Processes: Theory for Applications, Robert G. Gallager,
#' Cambridge University Press
#' 
#' @examples 
#' transMatr<-matrix(c(0.99,0.01,0.01,0.99),nrow=2,byrow=TRUE)
#' simpleMc<-new("markovchain", states=c("a","b"),
#'              transitionMatrix=transMatr)
#' expectedRewards(simpleMc,1,c(0,1))
#' @export
expectedRewards <- function(markovchain, n, rewards) {
  if (!is(markovchain, "markovchain")) {
    stop("markovchain must be a markovchain object")
  }
  nstates <- length(states(markovchain))
  if (length(n) != 1L || is.na(n) || !is.finite(n) ||
      n < 0 || n != floor(n)) {
    stop("n must be a non-negative integer")
  }
  if (!is.numeric(rewards) || length(rewards) != nstates || anyNA(rewards) ||
      any(!is.finite(rewards))) {
    stop("rewards must contain one finite numeric value for every state")
  }
  out <- .expectedRewardsRCpp(
    .rowStochasticMatrix(markovchain), as.integer(n), rewards)
  as.numeric(out)
}


#' Expected first passage Rewards for a set of states in a markovchain
#' 
#' @description Given a markovchain object and reward values for every state,
#' function calculates expected reward value for a set A of states after n 
#' steps. 
#'  
#' @usage expectedRewardsBeforeHittingA(markovchain, A, state, rewards, n)
#'  
#' @param markovchain the markovchain-class object
#' @param A set of states for first passage expected reward
#' @param state initial state
#' @param rewards vector depicting rewards coressponding to states
#' @param n no of steps of the process
#'  
#' @details The function returns the value of expected first passage 
#' rewards given rewards coressponding to every state, an initial state
#' and number of steps.
#'  
#' @return returns a expected reward (numerical value) as described above
#'  
#' @author Sai Bhargav Yalamanchi, Vandit Jain
#'  
#' @export
expectedRewardsBeforeHittingA <- function(markovchain, A, state, rewards, n) {
  if (!is(markovchain, "markovchain")) {
    stop("markovchain must be a markovchain object")
  }
  stateNames <- states(markovchain)
  if (!is.character(A) || length(A) < 1L || anyNA(A) ||
      any(!A %in% stateNames)) {
    stop("A must contain valid state names")
  }
  A <- unique(A)
  if (length(state) != 1L || is.na(state) || !state %in% stateNames) {
    stop("state must identify exactly one state")
  }
  if (state %in% A) {
    stop("the initial state must not belong to A")
  }
  if (!is.numeric(rewards) || length(rewards) != length(stateNames) ||
      anyNA(rewards) || any(!is.finite(rewards))) {
    stop("rewards must contain one finite numeric value for every state")
  }
  if (length(n) != 1L || is.na(n) || !is.finite(n) ||
      n < 0 || n != floor(n)) {
    stop("n must be a non-negative integer")
  }

  keep <- which(!stateNames %in% A)
  initial <- match(state, stateNames[keep])
  .expectedRewardsBeforeHittingARCpp(
    .rowStochasticMatrix(markovchain)[keep, keep, drop = FALSE],
    initial,
    rewards[keep],
    as.integer(n)
  )
}


#' Mean First Passage Time for irreducible Markov chains
#'
#' @description Given an irreducible (ergodic) markovchain object, this function
#'   calculates the expected number of steps to reach other states
#'
#' @param object the markovchain object
#' @param destination a character vector representing the states respect to
#'   which we want to compute the mean first passage time. Empty by default
#'
#' @details For an ergodic Markov chain it computes: 
#' \itemize{ 
#'   \item If destination is empty, the average first time (in steps) that takes
#'   the Markov chain to go from initial state i to j. (i, j) represents that 
#'   value in case the Markov chain is given row-wise, (j, i) in case it is given
#'   col-wise. 
#'   \item If destination is not empty, the average time it takes us from the 
#'   remaining states to reach the states in \code{destination} 
#' }
#'
#' @return a Matrix of the same size with the average first passage times if
#'   destination is empty, a vector if destination is not
#'
#' @author Toni Giorgino, Ignacio Cordón
#'
#' @references C. M. Grinstead and J. L. Snell. Introduction to Probability.
#' American Mathematical Soc., 2012.
#'
#' @examples
#' m <- matrix(1 / 10 * c(6,3,1,
#'                        2,3,5,
#'                        4,1,5), ncol = 3, byrow = TRUE)
#' mc <- new("markovchain", states = c("s","c","r"), transitionMatrix = m)
#' meanFirstPassageTime(mc, "r")
#'
#'
#' # Grinstead and Snell's "Oz weather" worked out example
#' mOz <- matrix(c(2,1,1,
#'                 2,0,2,
#'                 1,1,2)/4, ncol = 3, byrow = TRUE)
#'
#' mcOz <- new("markovchain", states = c("s", "c", "r"), transitionMatrix = mOz)
#' meanFirstPassageTime(mcOz)
#'
#' @export meanFirstPassageTime
setGeneric("meanFirstPassageTime", function(object, destination) {
  standardGeneric("meanFirstPassageTime")
})


setMethod("meanFirstPassageTime",  signature("markovchain", "missing"),
  function(object, destination) {
    destination = character()
    .meanFirstPassageTimeRcpp(object, destination)
  }
)

setMethod("meanFirstPassageTime",  signature("markovchain", "character"),
  function(object, destination) {
    states <- object@states
    incorrectStates <- setdiff(destination, states)
    
    if (length(incorrectStates) > 0)
      stop("Some of the states you provided in destination do not match states from the markovchain")

    result <- .meanFirstPassageTimeRcpp(object, destination)
    asVector <- as.vector(result)
    names(asVector) <- colnames(result)
    
    asVector
  }
)

#' Mean recurrence time
#'
#' @description Computes the expected time to return to a recurrent state
#'   in case the Markov chain starts there
#'
#' @usage meanRecurrenceTime(object)
#'
#' @param object the markovchain object
#'
#' @return For a Markov chain it outputs is a named vector with the expected 
#'   time to first return to a state when the chain starts there.
#'   States present in the vector are only the recurrent ones. If the matrix
#'   is ergodic (i.e. irreducible), then all states are present in the output
#'   and order is the same as states order for the Markov chain
#'
#' @author Ignacio Cordón
#'
#' @references C. M. Grinstead and J. L. Snell. Introduction to Probability.
#' American Mathematical Soc., 2012.
#'
#' @examples
#' m <- matrix(1 / 10 * c(6,3,1,
#'                        2,3,5,
#'                        4,1,5), ncol = 3, byrow = TRUE)
#' mc <- new("markovchain", states = c("s","c","r"), transitionMatrix = m)
#' meanRecurrenceTime(mc)
#'
#' @export meanRecurrenceTime
setGeneric("meanRecurrenceTime", function(object) {
  standardGeneric("meanRecurrenceTime")
})

setMethod("meanRecurrenceTime", "markovchain", function(object) {
  .meanRecurrenceTimeRcpp(object)
})


#' Mean absorption time
#'
#' @description Computes the expected number of steps to go from any of the
#'   transient states to any of the recurrent states. The Markov chain should
#'   have at least one transient state for this method to work
#'
#' @usage meanAbsorptionTime(object)
#'
#' @param object the markovchain object
#'
#' @return A named vector with the expected number of steps to go from a
#'   transient state to any of the recurrent ones
#'
#' @author Ignacio Cordón
#'
#' @references C. M. Grinstead and J. L. Snell. Introduction to Probability.
#' American Mathematical Soc., 2012.
#'
#' @examples
#' m <- matrix(c(1/2, 1/2, 0,
#'               1/2, 1/2, 0,
#'                 0, 1/2, 1/2), ncol = 3, byrow = TRUE)
#' mc <- new("markovchain", states = letters[1:3], transitionMatrix = m)
#' times <- meanAbsorptionTime(mc)
#'
#' @export meanAbsorptionTime
setGeneric("meanAbsorptionTime", function(object) {
  standardGeneric("meanAbsorptionTime")
})

setMethod("meanAbsorptionTime",  "markovchain", function(object) {
  .meanAbsorptionTimeRcpp(object)
})

#' Absorption probabilities
#'
#' @description Computes the absorption probability from each transient
#'   state to each recurrent one (i.e. the (i, j) entry or (j, i), in a 
#'   stochastic matrix by columns, represents the probability that the
#'   first not transient state we can go from the transient state i is j
#'   (and therefore we are going to be absorbed in the communicating
#'   recurrent class of j)
#'
#' @usage absorptionProbabilities(object)
#'
#' @param object the markovchain object
#'
#' @return A named vector with the expected number of steps to go from a
#'   transient state to any of the recurrent ones
#'
#' @author Ignacio Cordón
#'
#' @references C. M. Grinstead and J. L. Snell. Introduction to Probability.
#' American Mathematical Soc., 2012.
#'
#' @examples
#' m <- matrix(c(1/2, 1/2, 0,
#'               1/2, 1/2, 0,
#'                 0, 1/2, 1/2), ncol = 3, byrow = TRUE)
#' mc <- new("markovchain", states = letters[1:3], transitionMatrix = m)
#' absorptionProbabilities(mc)
#'
#' @export absorptionProbabilities
setGeneric("absorptionProbabilities", function(object) {
  standardGeneric("absorptionProbabilities")
})

setMethod("absorptionProbabilities",  "markovchain", function(object) {
  .absorptionProbabilitiesRcpp(object)
})


#' @title Check if a DTMC is regular
#' 
#' @description Function to check wether a DTCM is regular
# 
#' @details A Markov chain is regular if some of the powers of its matrix has all elements 
#'   strictly positive
#' 
#' @param object a markovchain object
#'
#' @return A boolean value
#'
#' @author Ignacio Cordón
#' @references Matrix Analysis. Roger A.Horn, Charles R.Johnson. 2nd edition. 
#'   Corollary 8.5.8, Theorem 8.5.9
#'
#' 
#' @examples 
#' P <- matrix(c(0.5,  0.25, 0.25,
#'               0.5,     0, 0.5,
#'               0.25, 0.25, 0.5), nrow = 3)
#' colnames(P) <- rownames(P) <- c("R","N","S")
#' ciao <- as(P, "markovchain")
#' is.regular(ciao)
#' 
#' @seealso \code{\link{is.irreducible}}
#' 
#' @exportMethod is.regular
setGeneric("is.regular", function(object) standardGeneric("is.regular"))

setMethod("is.regular", "markovchain", function(object) {
  .isRegularRcpp(object)
})


#' @title Check if a Markov chain is stochastically monotone
#' @description Verifies if the transition matrix of the Markov chain is stochastically monotone.
#' @param object A markovchain object or a transition matrix.
#' @return A boolean value.
#' @export
setGeneric("is.stochasticallyMonotone", function(object) standardGeneric("is.stochasticallyMonotone"))

#' @rdname is.stochasticallyMonotone
#' @aliases is.stochasticallyMonotone,markovchain-method
setMethod("is.stochasticallyMonotone", 
          signature(object = "markovchain"), 
          function(object) {
            return(.is_stochastically_monotone_cpp(.rowStochasticMatrix(object)))
          })

#' @rdname is.stochasticallyMonotone
#' @aliases is.stochasticallyMonotone,matrix-method
setMethod("is.stochasticallyMonotone", 
          signature(object = "matrix"), 
          function(object) {
            return(.is_stochastically_monotone_cpp(object))
          })
#' @rdname is.stochasticallyMonotone
#' @aliases is.stochasticallyMonotone,ANY-method
setMethod("is.stochasticallyMonotone", 
          signature(object = "ANY"), 
          function(object) {
            stop("must be a `markovchain` object or a matrix")
          })

#' Hitting probabilities for markovchain
#' 
#' @description Given a markovchain object,
#' this function calculates the probability of ever arriving from state i to j
#' 
#' @usage hittingProbabilities(object, targets = NULL,
#'   solver = c("direct", "bicgstab", "doubling"), tol = 1e-13, maxIter = 200)
#'
#' @param object the markovchain-class object
#' @param targets optional character vector of state names: only the hitting
#' probabilities \emph{towards} these states are computed. The default,
#' \code{NULL}, means all the states, which gives the full matrix as before.
#' Every target is handled independently of the others, so the work is
#' proportional to the number of targets: for a large chain, asking only for
#' the states of interest is much faster than computing the whole matrix and
#' subsetting it. Duplicated or unknown names are an error.
#' @param solver the method used to solve the linear system
#' \eqn{(I - Q) h = R} on the states whose probability is neither
#' structurally zero nor one (see Details). \code{"direct"}, the default, is
#' an LU factorisation: it costs \eqn{O(m^3)} once per target and is the most
#' accurate. \code{"bicgstab"} is an unpreconditioned BiCGSTAB iteration on
#' the sparse system: every iteration costs two sparse matrix-vector products
#' instead of a dense \eqn{O(m^3)} step, so it is the fastest choice on large
#' sparse chains, at the price of a looser residual. On a breakdown of the
#' iteration it restarts from the current residual, and if the breakdown
#' persists it switches to \code{"direct"} with a warning. \code{"doubling"} is the
#' doubled Neumann series used by versions up to 1.2, kept for reproducibility
#' of earlier results; it is also the automatic fallback of \code{"direct"} on
#' a numerically singular system.
#' @param tol relative residual at which the iterative solvers
#' (\code{"bicgstab"}, \code{"doubling"}) stop. Ignored by \code{"direct"},
#' except when it falls back to \code{"doubling"}.
#' @param maxIter maximum number of iterations of the iterative solvers: the
#' number of BiCGSTAB steps, or the number of squarings of the doubled Neumann
#' series. A warning is raised, and the current values returned, if the
#' requested \code{tol} is not reached within this many iterations.
#'
#' @details On each target the states are first split by graph reachability:
#' a state that cannot reach the target has probability zero, and one that can
#' reach the target but no closed class outside it has probability one. Only
#' the remaining states need the linear system that \code{solver} controls, so
#' on chains where that split already decides every state (an irreducible
#' chain, for instance) all three solvers do the same negligible amount of
#' work. The choice matters on chains with several closed classes, i.e. on
#' genuine absorption probabilities.
#'
#' @return a matrix of hitting probabilities. Entry \code{[i, j]} is the
#' probability of ever arriving from state \code{i} to state \code{j} (the
#' probability of returning, after at least one transition, on the diagonal);
#' for a chain with \code{byrow = FALSE} the matrix is transposed, as the
#' transition matrix is. With \code{targets}, only the columns (rows if
#' \code{byrow = FALSE}) of the targets are returned, in the order given, and
#' they coincide with those of the full matrix.
#'
#' @author Ignacio Cordón
#'
#' @references R. Vélez, T. Prieto, Procesos Estocásticos, Librería UNED, 2013
#'
#' H. A. van der Vorst (1992). Bi-CGSTAB: A Fast and Smoothly Converging
#' Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems.
#' \emph{SIAM Journal on Scientific and Statistical Computing}, 13(2), 631-644.
#'
#' @examples
#' M <- markovchain:::zeros(5)
#' M[1,1] <- M[5,5] <- 1
#' M[2,1] <- M[2,3] <- 1/2
#' M[3,2] <- M[3,4] <- 1/2
#' M[4,2] <- M[4,5] <- 1/2
#'
#' mc <- new("markovchain", transitionMatrix = M)
#' hittingProbabilities(mc)
#'
#' # only the probabilities of ever reaching the first state
#' hittingProbabilities(mc, targets = "1")
#'
#' # on a large sparse chain, the iterative solver avoids the dense products
#' hittingProbabilities(mc, targets = "1", solver = "bicgstab")
#'
#' @exportMethod hittingProbabilities
setGeneric("hittingProbabilities", function(object, targets = NULL,
                                            solver = c("direct", "bicgstab", "doubling"),
                                            tol = 1e-13, maxIter = 200)
  standardGeneric("hittingProbabilities"))

setMethod("hittingProbabilities", "markovchain", function(object, targets = NULL,
                                                          solver = c("direct", "bicgstab", "doubling"),
                                                          tol = 1e-13, maxIter = 200) {
  allStates <- object@states
  if (is.null(targets)) {
    idx <- seq_along(allStates)
  } else {
    if (!is.character(targets) || length(targets) < 1L || anyNA(targets))
      stop("targets must be NULL or a non-empty character vector of state names with no missing values.")
    if (anyDuplicated(targets))
      stop("targets must not contain duplicate state names.")
    idx <- match(targets, allStates)
    if (anyNA(idx))
      stop("Unknown state(s) in targets: ", paste(targets[is.na(idx)], collapse = ", "))
  }

  solver <- match.arg(solver)
  # Keep in sync with the HittingSolver enum in src/probabilistic.cpp.
  solverCode <- switch(solver, direct = 0L, bicgstab = 1L, doubling = 2L)

  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0)
    stop("tol must be a single positive number.")
  if (!is.numeric(maxIter) || length(maxIter) != 1L || is.na(maxIter) ||
      maxIter < 1 || maxIter != floor(maxIter))
    stop("maxIter must be a single positive integer.")

  .hittingProbabilitiesRcpp(object, as.integer(idx), solverCode,
                            as.numeric(tol), as.integer(maxIter))
})



#' Mean num of visits for markovchain, starting at each state
#' 
#' @description Given a markovchain object, this function calculates 
#' a matrix where the element (i, j) represents the expect number of visits
#' to the state j if the chain starts at i (in a Markov chain by columns it
#' would be the element (j, i) instead)
#' 
#' @usage meanNumVisits(object)
#' 
#' @param object the markovchain-class object
#' 
#' @return a matrix with the expect number of visits to each state
#' 
#' @author Ignacio Cordón
#' 
#' @references R. Vélez, T. Prieto, Procesos Estocásticos, Librería UNED, 2013
#' 
#' @examples
#' M <- markovchain:::zeros(5)
#' M[1,1] <- M[5,5] <- 1
#' M[2,1] <- M[2,3] <- 1/2
#' M[3,2] <- M[3,4] <- 1/2
#' M[4,2] <- M[4,5] <- 1/2
#' 
#' mc <- new("markovchain", transitionMatrix = M)
#' meanNumVisits(mc)
#' 
#' @exportMethod meanNumVisits
setGeneric("meanNumVisits", function(object) standardGeneric("meanNumVisits"))

setMethod("meanNumVisits", "markovchain", function(object) {
  .minNumVisitsRcpp(object)
})


setMethod(
  "steadyStates",
  "markovchain", 
  function(object) {
    .steadyStatesRcpp(object)
  }
)


# Internal helper: the properties shown by summary(object, details = TRUE),
# in the spirit of the printout of a PyDTMC MarkovChain. Quantities that
# need a unique stationary distribution (or an irreducible chain) are NA
# when they are not defined.
.summaryDetails <- function(object) {
  P <- .rowStochasticMatrix(object)
  n <- nrow(P)
  tol <- sqrt(.Machine$double.eps)
  safe <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) NA)

  classes <- communicatingClasses(object)
  recurrent <- recurrentClasses(object)
  irreducible <- is.irreducible(object)
  absorbing <- absorbingStates(object)
  per <- if (irreducible) safe(period(object)) else NA_integer_
  uniqueStationary <- length(recurrent) == 1L

  list(
    size = n,
    rank = qr(P)$rank,
    classes = length(classes),
    recurrentClasses = length(recurrent),
    transientClasses = length(classes) - length(recurrent),
    irreducible = irreducible,
    period = per,
    regular = isTRUE(safe(is.regular(object))),
    absorbingChain = length(absorbing) > 0L && all(lengths(recurrent) == 1L),
    reversible = if (irreducible) isTRUE(safe(is.reversible(object))) else NA,
    stochasticallyMonotone = isTRUE(safe(is.stochasticallyMonotone(object))),
    symmetric = isTRUE(all.equal(P, t(P), tolerance = tol, check.attributes = FALSE)),
    entropyRate = if (uniqueStationary) safe(entropyRate(object)) else NA_real_,
    slem = if (irreducible) safe(slem(object)) else NA_real_,
    spectralGap = if (irreducible) safe(spectralGap(object)) else NA_real_,
    kemenyConstant = if (irreducible) safe(kemenyConstant(object)) else NA_real_
  )
}

.printSummaryDetails <- function(d) {
  yn <- function(x) if (is.na(x)) "not defined" else if (x) "yes" else "no"
  num <- function(x) if (is.na(x)) "not defined" else format(signif(x, 6))
  rows <- c(
    "Size" = as.character(d$size),
    "Rank" = as.character(d$rank),
    "Communicating classes" = paste0(d$classes, " (", d$recurrentClasses,
                                     " recurrent, ", d$transientClasses,
                                     " transient)"),
    "Irreducible" = yn(d$irreducible),
    "Period" = if (is.na(d$period)) "not defined" else as.character(d$period),
    "Regular (ergodic)" = yn(d$regular),
    "Absorbing chain" = yn(d$absorbingChain),
    "Reversible" = yn(d$reversible),
    "Stochastically monotone" = yn(d$stochasticallyMonotone),
    "Symmetric" = yn(d$symmetric),
    "Entropy rate (bits)" = num(d$entropyRate),
    "SLEM" = num(d$slem),
    "Spectral gap" = num(d$spectralGap),
    "Kemeny constant" = num(d$kemenyConstant)
  )
  cat("Further properties:", "\n")
  w <- max(nchar(names(rows)))
  for (k in seq_along(rows)) {
    cat(" ", formatC(names(rows)[k], width = -w), ":", rows[[k]], "\n")
  }
}

#' @exportMethod summary
setGeneric("summary")

# summary method for markovchain class
# lists: closed, transient classes, irreducibility, absorbint, transient states
setMethod("summary", signature(object = "markovchain"),
  function(object, details = FALSE, ...){
    if (length(details) != 1L || !is.logical(details) || is.na(details)) {
      stop("details must be TRUE or FALSE.")
    }
    
    # list of closed, recurrent and transient classes
    outs <- .summaryKernelRcpp(object)
    
    # display name of the markovchain object
    cat(object@name," Markov chain that is composed by:", "\n")
    
    # number of closed classes
    check <- length(outs$closedClasses)
    
    cat("Closed classes:","\n")
    
    # display closed classes
    if(check == 0) cat("NONE", "\n") else {
      for(i in 1:check) cat(outs$closedClasses[[i]], "\n")
    }
    
    # number of recurrent classes
    check <- length(outs$recurrentClasses)
    
    cat("Recurrent classes:", "\n")
    
    # display recurrent classes
    if(check == 0) cat("NONE", "\n") else {
      cat("{")
      cat(outs$recurrentClasses[[1]], sep = ",")
      cat("}")
      if(check > 1) {
        for(i in 2:check) {
          cat(",{")
          cat(outs$recurrentClasses[[i]], sep = ",")
          cat("}")
        }
      }
      cat("\n")
    }
    
    # number of transient classes
    check <- length(outs$transientClasses)
    
    cat("Transient classes:","\n")
    
    # display transient classes
    if(check == 0) cat("NONE", "\n") else {
      cat("{")
      cat(outs$transientClasses[[1]], sep = ",")
      cat("}")
      if(check > 1) { 
        for(i in 2:check) {
          cat(",{")
          cat(outs$transientClasses[[i]], sep = ",")
          cat("}")
        }
      }
      cat("\n")
    }
    
    # bool to say about irreducibility of markovchain
    irreducibility <- is.irreducible(object)
    
    if(irreducibility) 
      cat("The Markov chain is irreducible", "\n") 
    else cat("The Markov chain is not irreducible", "\n")
    
    # display absorbing states
    check <- absorbingStates(object)
    if(length(check) == 0) check <- "NONE"
    cat("The absorbing states are:", check )
    cat("\n")
    
    # optional block of further properties, printed after the classic output
    # so that summary(object) itself is unchanged
    if (details) {
      outs$details <- .summaryDetails(object)
      .printSummaryDetails(outs$details)
    }
    
    # return outs
    # useful when user will assign the value returned
    invisible(outs) 
  }
)

# Validate and convert a partition to C++ indices.
#
# This internal helper checks that `partition` is a named, exhaustive and
# mutually exclusive partition of the state space, then converts state names to
# zero-based integer indices for the C++ backend.
.get_partition_indices <- function(state_names, partition) {
  if (!is.list(partition) || length(partition) < 1L) {
    stop("Invalid partition: partition must be a non-empty list.")
  }
  if (is.null(names(partition)) || any(!nzchar(names(partition)))) {
    stop("Invalid partition: partition must be a named list.")
  }
  if (anyDuplicated(names(partition))) {
    stop("Invalid partition: macro-state names must be unique.")
  }

  part_idx <- lapply(partition, function(x) {
    if (!is.character(x) || length(x) < 1L) {
      stop("Invalid partition: each macro-state must contain at least one state name.")
    }
    idx <- match(x, state_names)
    if (any(is.na(idx))) {
      stop("Invalid partition: Some states in the partition do not exist in the Markov chain.")
    }
    as.integer(idx - 1L)
  })

  all_idx <- unlist(part_idx, use.names = FALSE)
  if (length(all_idx) != length(state_names)) {
    stop("Invalid partition: The partition must contain all states of the Markov chain exactly once (no duplicates, no omissions).")
  }
  if (length(unique(all_idx)) != length(state_names)) {
    stop("Invalid partition: The partition must contain all states of the Markov chain exactly once (no duplicates, no omissions).")
  }

  part_idx
}

#' Check exact lumpability of a Markov chain
#'
#' @description Verifies the strong lumpability condition with respect to a
#' partition of the state space. For every pair of macro-states, all micro-states
#' in the same source macro-state must have the same total probability of moving
#' to the destination macro-state.
#'
#' @param object A \code{markovchain} object.
#' @param partition A named list of character vectors defining macro-states.
#' @param tol Non-negative numerical tolerance for equality checks.
#' @return A logical value.
#' @references Kemeny, J. G. and Snell, J. L. (1960). \emph{Finite Markov Chains}.
#' @export
setGeneric("is.lumpable", function(object, partition, tol = 1e-10) standardGeneric("is.lumpable"))

#' @rdname is.lumpable
#' @aliases is.lumpable,markovchain-method
setMethod("is.lumpable", signature(object = "markovchain"),
          function(object, partition, tol = 1e-10) {
            part_idx <- .get_partition_indices(states(object), partition)
            P <- object@transitionMatrix
            # The C++ backend checks row-wise transition probabilities.  If the
            # object stores probabilities by column, transpose the matrix so the
            # lumpability condition is still evaluated on outgoing probabilities.
            if (!object@byrow) {
              P <- t(P)
            }
            .is_lumpable_cpp(P, part_idx, tol)
          })

#' Aggregate a Markov chain over a partition
#'
#' @description Coarsens a Markov chain to a reduced state space. By default the
#' function requires exact lumpability. With \code{force = TRUE}, it performs an
#' approximate aggregation using stationary weights when available and arithmetic
#' averages for macro-states with zero stationary mass.
#'
#' @param object A \code{markovchain} object.
#' @param partition A named list of character vectors defining macro-states.
#' @param force If \code{FALSE}, stop unless the chain is exactly lumpable. If
#' \code{TRUE}, return a weighted approximate lumping.
#' @return A \code{markovchain} object on the macro-state space.
#' @export
setGeneric("lump", function(object, partition, force = FALSE) standardGeneric("lump"))

#' @rdname lump
#' @aliases lump,markovchain-method
setMethod("lump", signature(object = "markovchain"),
          function(object, partition, force = FALSE) {
            part_idx <- .get_partition_indices(states(object), partition)

            P <- object@transitionMatrix
            # Work internally with the usual row-stochastic convention.  The
            # returned object is also row-stochastic, independently of the input
            # storage orientation.
            if (!object@byrow) {
              P <- t(P)
            }

            if (!force && !.is_lumpable_cpp(P, part_idx, 1e-10)) {
              stop("The Markov chain is not exactly lumpable. Use force = TRUE to perform an approximate weighted lumping.")
            }

            st <- steadyStates(object)
            # steadyStates() returns one distribution per row for row-stored
            # chains and one per column otherwise; put them in rows.
            if (!object@byrow) st <- t(st)
            if (nrow(st) > 0L) {
              # If several stationary distributions are returned, average them
              # to obtain deterministic non-negative aggregation weights.
              w <- colMeans(st)
            } else {
              w <- rep(1 / ncol(object@transitionMatrix), ncol(object@transitionMatrix))
            }

            P_lumped <- .lump_cpp(P, part_idx, as.numeric(w))
            dimnames(P_lumped) <- list(names(partition), names(partition))

            new("markovchain",
                states = names(partition),
                transitionMatrix = P_lumped,
                byrow = TRUE,
                name = paste(object@name, "(Lumped)"))
          })

#' Automatically aggregate a Markov chain by spectral clustering
#'
#' @description Finds an approximate partition by clustering the leading right
#' eigenvectors of the transition matrix, then returns the forced lumping over
#' that partition. This is a heuristic for approximate lumping/metastable
#' aggregation, not a proof of exact lumpability.
#'
#' @param object A \code{markovchain} object.
#' @param k Number of macro-states to discover.
#' @return A list with \code{partition} and \code{lumped_chain}.
#' @export
setGeneric("autoLump", function(object, k) standardGeneric("autoLump"))

#' @rdname autoLump
#' @aliases autoLump,markovchain-method
setMethod("autoLump", signature(object = "markovchain"),
          function(object, k) {
            # eigenvectors must be those of the row-stochastic matrix, whatever
            # the storage orientation of the chain
            P <- .rowStochasticMatrix(object)
            state_names <- states(object)
            n <- nrow(P)

            if (length(k) != 1L || is.na(k) || k != as.integer(k)) {
              stop("k must be a single integer.")
            }
            k <- as.integer(k)
            if (k <= 1L || k >= n) {
              stop("The number of macro-states 'k' must be between 2 and the number of states - 1.")
            }

            eig <- eigen(P)
            ord <- order(Mod(eig$values), decreasing = TRUE)
            V_eig <- Re(eig$vectors[, ord[seq_len(k)], drop = FALSE])

            # Make the example deterministic without permanently changing the
            # user's random-number stream.
            old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
              get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
            } else {
              NULL
            }
            on.exit({
              if (is.null(old_seed)) {
                if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                  rm(".Random.seed", envir = .GlobalEnv)
                }
              } else {
                assign(".Random.seed", old_seed, envir = .GlobalEnv)
              }
            }, add = TRUE)
            set.seed(42)
            clust <- stats::kmeans(V_eig, centers = k, nstart = 10)

            partition <- stats::setNames(vector("list", k), paste0("Macro_", seq_len(k)))
            for (i in seq_len(k)) {
              partition[[i]] <- state_names[clust$cluster == i]
            }

            list(
              partition = partition,
              lumped_chain = lump(object, partition, force = TRUE)
            )
          })
