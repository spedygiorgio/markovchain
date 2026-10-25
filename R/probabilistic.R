# Overview: most functions in this file are thin R wrappers (generic + method) around the C++ code in
# src/probabilistic.cpp. The R layer gives them a public name, validates the arguments and puts the
# matrix in the row-stochastic orientation; the algorithms (graph search, linear systems, ...) are in C++.
# Legend for the vocabulary used here: a 'class' is a set of states that communicate (can reach each
# other); a class is 'recurrent' if it cannot be left, otherwise 'transient'; an 'absorbing' state
# is a recurrent class made of one state.
# given a markovchain object is it possible to reach goal state from 
# a given state

#' @name is.accessible
#' @title Check whether a state can be reached from another
#' @description Checks whether state `to` is accessible from state `from`,
#'   i.e. whether the chain started in `from` reaches `to` with positive
#'   probability in some number of steps.
#'
#' @param object A [`markovchain`][markovchain-class] object.
#' @param from Name of the starting state.
#' @param to Name of the target state.
#' 
#' @return `TRUE` if `to` is accessible from `from`, `FALSE` otherwise.
#' 
#' @references James Montgomery, University of Madison
#' 
#' @author Giorgio Spedicato, Ignacio Cordón
#' @seealso [is.irreducible()], [communicatingClasses()]
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
# S4 generic: `is.accessible(chain, from, to)`. S4 methods are chosen according to the CLASS of ALL
# the arguments listed in the signature.
setGeneric("is.accessible", function(object, from, to) standardGeneric("is.accessible"))

# Method 1: two state names -> TRUE/FALSE, 'is `to` reachable from `from`'.
setMethod("is.accessible", c("markovchain", "character", "character"), 
  function(object, from, to) {
    # O(n²) procedure to see if to state is reachable starting at from state
    return(.isAccessibleRcpp(object, from, to))
  }
)

# Method 2: no state names -> the whole reachability matrix (see zzz_reachability.R, which replaces
# this method with a version that handles byrow = FALSE).
setMethod("is.accessible", c("markovchain", "missing", "missing"), 
  function(object, from, to) {
    .reachabilityMatrixRcpp(object)
  }
)

# a markov chain is irreducible if it is composed of only one communicating class

#' @name is.irreducible
#' @title Check whether a Markov chain is irreducible
#' @description A Markov chain is irreducible if all its states communicate,
#'   i.e. form a single communicating class.
#' @param object A [`markovchain`][markovchain-class] object.
#' 
#' @details Irreducibility does not imply aperiodicity: an irreducible chain
#'   on finitely many states has a unique stationary distribution, but it
#'   converges to it from every initial state only if it is also aperiodic
#'   (see [period()] and [is.regular()]).
#' @return `TRUE` if the chain is irreducible, `FALSE` otherwise.
#' 
#' @references Feres, Matlab listings for Markov Chains.
#' @author Giorgio Spedicato
#' 
#' @seealso [communicatingClasses()], [is.regular()], [period()]
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
# Irreducible: the whole chain is a single communicating class (every state reaches every other).
setGeneric("is.irreducible", function(object) standardGeneric("is.irreducible"))

setMethod("is.irreducible", "markovchain", function(object) {
  .isIrreducibleRcpp(object)
})


# what this function will do?
# It calculates the probability to go from given state
# to all other states in k steps
# k varies from 1 to n

#' @name firstPassage
#' @title First passage time distribution
#' @description Distribution of the first passage times of a Markov chain
#'   started in `state`: for every state \eqn{j} and step \eqn{k}, the
#'   probability that the chain enters \eqn{j} for the first time at step
#'   \eqn{k}.
#'
#' @param object A [`markovchain`][markovchain-class] object.
#' @param state Name of the initial state.
#' @param n Number of steps, a positive integer.
#'
#' @details For \eqn{j} equal to `state` the probabilities are those of the
#'   first return. The algorithm follows the Matlab listings of Feres.
#' @return A matrix with `n` rows (steps \eqn{1, \dots, n}) and one column per
#'   state: entry \eqn{(k, j)} is the probability that the first passage to
#'   \eqn{j} happens at step \eqn{k}. The column sums are at most one; the
#'   remainder is the probability of a first passage after step `n`, or never.
#'
#' @references Feres, R. Notes for Math 450: Matlab listings for Markov chains.
#'   Washington University in St. Louis.
#'
#' @author Giorgio Spedicato
#' @family first passage and hitting
#'
#' @examples
#' simpleMc <- new("markovchain", states = c("a", "b"),
#'                 transitionMatrix = matrix(c(0.4, 0.6, 0.3, 0.7),
#'                                           nrow = 2, byrow = TRUE))
#' firstPassage(simpleMc, "b", 20)
#' @export
# firstPassage(): probability that the chain started in `state` reaches each state for the FIRST time
# exactly at step 1, 2, ..., n. Row k of the result is the distribution of the first passage at time k.
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

  # the C++ code works on the row-stochastic matrix and returns an n x (number of states) matrix
  outMatr <- .firstpassageKernelRcpp(
    P = .rowStochasticMatrix(object),
    i = match(state, stateNames),
    n = as.integer(n)
  )
  colnames(outMatr) <- stateNames
  rownames(outMatr) <- seq_len(nrow(outMatr))
  outMatr
}


#' First passage time distribution of a set of states
#'
#' @description Distribution of the first passage time into a set of states:
#'   the probability that a Markov chain started in `state` enters `set` for
#'   the first time at step \eqn{k = 1, \dots, n}.
#'
#' @param object A [`markovchain`][markovchain-class] object.
#' @param state Name of the initial state.
#' @param set Names of the states of the target set.
#' @param n Number of steps, a positive integer.
#'
#' @return A one-column matrix (column `"set"`) with `n` rows: row \eqn{k}
#'   is the probability that the first passage into `set` happens at step
#'   \eqn{k}.
#' @references
#' Feres, R. Notes for Math 450: Matlab listings for Markov chains.
#' Washington University in St. Louis.
#'
#' Gallager, R. G. 6.262 Discrete Stochastic Processes, course notes,
#' chapter 5. MIT OpenCourseWare.
#'
#' @author Vandit Jain
#' @family first passage and hitting
#' @examples
#' statesNames <- c("a", "b", "c")
#' markovB <- new("markovchain", states = statesNames, transitionMatrix =
#'   matrix(c(0.2, 0.5, 0.3,
#'            0, 1, 0,
#'            0.1, 0.8, 0.1), nrow = 3, byrow = TRUE,
#'          dimnames = list(statesNames, statesNames)))
#' firstPassageMultiple(markovB, "a", c("b", "c"), 4)
#' @export
# firstPassageMultiple(): same idea, but the target is a SET of states: probability of entering the
# set for the first time at step k.
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

  # `set` is a single absorbing target for the C++ code: the states of the set are
  # identified by position (duplicates removed)
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
#' @title Structural analysis of a Markov chain
#' @description Classification of the states of a Markov chain from the graph
#'   of its possible transitions: communicating, recurrent and transient
#'   classes and states, absorbing states, the canonical form of the
#'   transition matrix, and the period.
#' 
#' @param object A \code{markovchain} object.
#' 
#' @return
#' \describe{
#'   \item{\code{period}}{returns an integer, the period of the Markov chain
#'     (the greatest common divisor of the lengths of its cycles); it is
#'     defined for an irreducible chain (otherwise a warning is given and 0
#'     is returned)}
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
#' transientClasses(mc)
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
# Communicating classes, recurrent classes, transient classes and states, absorbing states and
# canonical form: each is computed in C++ from the graph of possible transitions.
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


#' Committor probabilities
#'
#' @description Probability that a Markov chain hits the set of states `A`
#'   before the set `B` (the committor function), from every initial state or
#'   from one.
#'
#' @param object A [`markovchain`][markovchain-class] object.
#' @param A,B Disjoint sets of states, given as indices (positions in
#'   `states(object)`).
#' @param p Optional index of the initial state. If omitted, the probabilities
#'   for all initial states are returned.
#'
#' @details The committor \eqn{q} solves the linear system
#'   \eqn{q_i = \sum_j P_{ij} q_j} for the states outside \eqn{A \cup B},
#'   with \eqn{q = 1} on `A` and \eqn{q = 0} on `B`. The system has a unique
#'   solution when `A` or `B` can be reached from every state.
#'
#' @return A numeric vector with the probability for each initial state, or
#'   a single number if `p` is given.
#'
#' @family first passage and hitting
#' @examples
#' transMatr <- matrix(c(0, 0, 0, 1, 0.5,
#'                       0.5, 0, 0, 0, 0,
#'                       0.5, 0, 0, 0, 0,
#'                       0, 0.2, 0.4, 0, 0,
#'                       0, 0.8, 0.6, 0, 0.5), nrow = 5)
#' object <- new("markovchain", states = c("a", "b", "c", "d", "e"),
#'               transitionMatrix = transMatr)
#' # probability of reaching "e" before "c", from each state
#' committorAB(object, 5, 3)
#' committorAB(object, 5, 3, p = 1)
#' @export
# committorAB(): probability of hitting the set A before the set B (committor probability).
# It solves the linear system  (P - I) q = 0  for the states outside A and B, with q = 1 on one set
# and 0 on the other; the code below rewrites the rows of A and B to impose these boundary conditions.
committorAB <- function(object, A, B, p = 1) {
  if (!is(object, "markovchain")) {
    stop("object must be a markovchain object")
  }

  nstates <- length(object@states)
  # A, B are given as positions of states (1-based) in the order of states(object)
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

  # without `p` the committor of every state is returned, otherwise only that of state p
  return_all <- missing(p)
  if (!return_all &&
      (length(p) != 1L || is.na(p) || !is.finite(p) ||
       p != floor(p) || p < 1L || p > nstates)) {
    stop("please provide a valid initial state")
  }

  # coefficient matrix of the system, with the boundary conditions q[A] = 1 and q[B] = 0.
  coefficient <- .rowStochasticMatrix(object) - diag(nstates)
  coefficient[A, ] <- 0
  coefficient[cbind(A, A)] <- 1       # row of a state of A: the equation q_i = 1
  coefficient[B, ] <- 0
  coefficient[cbind(B, B)] <- 1       # row of a state of B: the equation q_i = 0
  rhs <- numeric(nstates)             # 0 outside A, and for B
  rhs[A] <- 1
  out <- solve(coefficient, rhs)      # q_i = sum_j P_ij q_j for every other state

  if (return_all) out else out[as.integer(p)]
}


#' Expected cumulative reward of a Markov chain
#'
#' @description Expected total reward collected by a Markov chain over the
#'   times \eqn{0, 1, \dots, n}, for every initial state, when state \eqn{i}
#'   pays `rewards[i]` at each visit.
#'
#' @param markovchain A [`markovchain`][markovchain-class] object.
#' @param n Number of steps, a non-negative integer.
#' @param rewards Numeric vector with the reward of each state, in the order
#'   of `states(markovchain)`.
#'
#' @details The values solve the dynamic programming recursion
#'   \eqn{v_0 = r}, \eqn{v_k = r + P v_{k-1}} (Gallager 2013, Section 3.5),
#'   so that \eqn{v_n = \sum_{k=0}^{n} P^k r}.
#'
#' @return A numeric vector with one value per initial state.
#'
#' @author Vandit Jain
#'
#' @references Gallager, R. G. (2013). Stochastic Processes: Theory for
#'   Applications. Cambridge University Press.
#' @family rewards
#' @examples
#' transMatr <- matrix(c(0.99, 0.01, 0.01, 0.99), nrow = 2, byrow = TRUE)
#' simpleMc <- new("markovchain", states = c("a", "b"),
#'                 transitionMatrix = transMatr)
#' expectedRewards(simpleMc, 1, c(0, 1))
#' @export
# expectedRewards(): each state i gives a reward r_i each time it is visited; the function returns,
# for each starting state, the expected total reward after n steps = sum_{k<=n} P^k r.
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


#' Expected reward collected before hitting a set of states
#'
#' @description Expected total reward collected by a Markov chain started in
#'   `state` at the steps \eqn{1, \dots, n} it spends outside the set `A`
#'   before entering it, when state \eqn{i} pays `rewards[i]` at each visit.
#'
#' @param markovchain A [`markovchain`][markovchain-class] object.
#' @param A Names of the states of the target set.
#' @param state Name of the initial state, not in `A`.
#' @param rewards Numeric vector with the reward of each state, in the order
#'   of `states(markovchain)` (the rewards of the states in `A` are not used).
#' @param n Number of steps, a non-negative integer.
#'
#' @details With \eqn{Q} the transition matrix restricted to the states
#'   outside `A`, the value is \eqn{\sum_{k=1}^{n} (Q^k r)_{s}}, with
#'   \eqn{s} the initial state: the reward of the initial state at time 0 is
#'   not included.
#'
#' @return A number, the expected reward.
#'
#' @author Sai Bhargav Yalamanchi, Vandit Jain
#' @family rewards
#' @examples
#' transMatr <- matrix(c(0.5, 0.3, 0.2,
#'                       0.2, 0.5, 0.3,
#'                       0, 0, 1), nrow = 3, byrow = TRUE)
#' mc <- new("markovchain", states = c("a", "b", "c"),
#'           transitionMatrix = transMatr)
#' # reward 1 per step spent in "a" or "b" before absorption in "c"
#' expectedRewardsBeforeHittingA(mc, A = "c", state = "a",
#'                               rewards = c(1, 1, 0), n = 100)
#' @export
# expectedRewardsBeforeHittingA(): same as above but the chain is stopped when it enters the set A.
# The states in A are removed from the matrix, so only the transitions among the remaining states count.
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

  # Keep only the states outside A: the reward is collected only while the chain has not
  # entered A, and `initial` is the position of the starting state among the kept ones.
  keep <- which(!stateNames %in% A)
  initial <- match(state, stateNames[keep])
  .expectedRewardsBeforeHittingARCpp(
    .rowStochasticMatrix(markovchain)[keep, keep, drop = FALSE],
    initial,
    rewards[keep],
    as.integer(n)
  )
}


#' Mean first passage times of an irreducible Markov chain
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
#' @family first passage and hitting
#' @export meanFirstPassageTime
# Mean first passage time: expected number of steps to go from i to a destination state.
setGeneric("meanFirstPassageTime", function(object, destination) {
  standardGeneric("meanFirstPassageTime")
})


# No destination: matrix of the mean first passage times between every pair of states
setMethod("meanFirstPassageTime",  signature("markovchain", "missing"),
  function(object, destination) {
    destination = character()
    .meanFirstPassageTimeRcpp(object, destination)
  }
)

# With a destination (a set of states): vector with the mean time from each other state
# to reach the set
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
#' @family first passage and hitting
#' @export meanRecurrenceTime
# Mean recurrence time of a state: expected number of steps to return to it (1 / stationary prob.).
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
#' @family first passage and hitting
#' @export meanAbsorptionTime
# Mean absorption time: expected number of steps before being absorbed, from each transient state.
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
#' @family first passage and hitting
#' @export absorptionProbabilities
# Absorption probabilities: probability of ending in each absorbing state, from each transient state.
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
# Regular: some power of the matrix has all strictly positive entries.
setGeneric("is.regular", function(object) standardGeneric("is.regular"))

setMethod("is.regular", "markovchain", function(object) {
  .isRegularRcpp(object)
})


#' @title Check whether a Markov chain is stochastically monotone
#' @description A Markov chain on ordered states is stochastically monotone
#'   if starting from a higher state never makes the next state
#'   stochastically smaller: for every pair of consecutive states \eqn{i},
#'   \eqn{i + 1} and every \eqn{j}, \eqn{\sum_{l \le j} P_{i+1, l} \le
#'   \sum_{l \le j} P_{i, l}}.
#' @param object A [`markovchain`][markovchain-class] object, or a
#'   row-stochastic transition matrix.
#' @details The states are taken in the order of the transition matrix, and
#'   the inequalities are checked up to a tolerance of \eqn{10^{-10}}.
#' @return `TRUE` if the chain is stochastically monotone, `FALSE` otherwise.
#' @references Keilson, J. and Kester, A. (1977). Monotone matrices and
#'   monotone Markov processes. Stochastic Processes and their Applications
#'   5(3), 231-241.
#' @examples
#' # a birth-death chain is monotone
#' P <- matrix(c(0.7, 0.3, 0,
#'               0.2, 0.5, 0.3,
#'               0, 0.4, 0.6), nrow = 3, byrow = TRUE,
#'             dimnames = list(c("low", "mid", "high"), c("low", "mid", "high")))
#' is.stochasticallyMonotone(new("markovchain", transitionMatrix = P))
#' # reversing the order of the rows breaks monotonicity
#' is.stochasticallyMonotone(P[3:1, ])
#' @export
# Stochastically monotone: higher current states imply stochastically higher next states.
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

#' Hitting probabilities of a Markov chain
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
#' @family first passage and hitting
#' @exportMethod hittingProbabilities
# hittingProbabilities(): probability that the chain started at state i ever reaches the target state.
# Solvers: 'direct' (LU factorisation), 'bicgstab' (iterative, for large sparse chains) and
# 'doubling' (older method); the numeric code is passed to C++.
setGeneric("hittingProbabilities", function(object, targets = NULL,
                                            solver = c("direct", "bicgstab", "doubling"),
                                            tol = 1e-13, maxIter = 200)
  standardGeneric("hittingProbabilities"))

setMethod("hittingProbabilities", "markovchain", function(object, targets = NULL,
                                                          solver = c("direct", "bicgstab", "doubling"),
                                                          tol = 1e-13, maxIter = 200) {
  # Step 1: choose the target states (all of them by default) and validate.
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

  # Step 2: choose and validate the solver, then call C++.
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



#' Mean number of visits to each state
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
#' @family first passage and hitting
#' @exportMethod meanNumVisits
# Mean number of visits to each transient state before leaving the transient states.
setGeneric("meanNumVisits", function(object) standardGeneric("meanNumVisits"))

setMethod("meanNumVisits", "markovchain", function(object) {
  .minNumVisitsRcpp(object)
})


# The steady states (stationary distributions, one for each recurrent class) are computed in C++.
setMethod(
  "steadyStates",
  "markovchain", 
  function(object) {
    .steadyStatesRcpp(object)
  }
)


# Builds the extra block of properties shown by summary(details = TRUE). `safe()` runs an expression
# and returns NA instead of failing, because some properties are not defined for every chain.
# Internal helper: the properties shown by summary(object, details = TRUE),
# in the spirit of the printout of a PyDTMC MarkovChain. Quantities that
# need a unique stationary distribution (or an irreducible chain) are NA
# when they are not defined.
.summaryDetails <- function(object) {
  P <- .rowStochasticMatrix(object)
  n <- nrow(P)
  tol <- sqrt(.Machine$double.eps)
  # safe(): evaluates an expression and returns NA if it fails
  safe <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) NA)

  classes <- communicatingClasses(object)
  recurrent <- recurrentClasses(object)
  irreducible <- is.irreducible(object)
  absorbing <- absorbingStates(object)
  per <- if (irreducible) safe(period(object)) else NA_integer_
  # exactly one recurrent class <=> unique stationary distribution (needed for the entropy rate)
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

# Prints the block built by .summaryDetails() as two aligned columns.
.printSummaryDetails <- function(d) {
  # yn(), num(): format a logical / numeric property, "not defined" for NA
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

# summary() for chains: prints closed classes, recurrent classes, transient classes, irreducibility and
# absorbing states. The classic output is unchanged; details = TRUE adds more properties.
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

# Lumping: merging several states into one macro-state. A partition is a named list: each element is a
# vector of state names. This helper checks that the partition covers every state exactly once, then
# converts names into 0-based positions (C++ counts from 0, R from 1).
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

  # each block -> its states as 0-based positions in state_names (what the C++ code expects)
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

  # exhaustive and mutually exclusive: every state appears exactly once across the blocks
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
#' @description Checks whether a Markov chain is (strongly) lumpable with
#'   respect to a partition of its states: for every pair of blocks
#'   \eqn{A, B} of the partition, all the states of \eqn{A} have the same
#'   probability of moving into \eqn{B} (Kemeny and Snell 1960, Theorem
#'   6.3.2). The aggregated process is then a Markov chain, for every initial
#'   distribution.
#'
#' @param object A [`markovchain`][markovchain-class] object.
#' @param partition A named list of character vectors, the blocks of the
#'   partition (macro-states); every state must belong to exactly one block.
#' @param tol Non-negative tolerance for the equality of the probabilities
#'   (default `1e-10`).
#' @return `TRUE` if the chain is lumpable with respect to `partition`,
#'   `FALSE` otherwise.
#' @references Kemeny, J. G. and Snell, J. L. (1960). Finite Markov Chains.
#'   Van Nostrand.
#' @family state aggregation
#' @examples
#' statesNames <- c("a", "b", "c", "d")
#' P <- matrix(c(0.1, 0.2, 0.3, 0.4,
#'               0.2, 0.1, 0.4, 0.3,
#'               0.3, 0.3, 0.2, 0.2,
#'               0.2, 0.4, 0.1, 0.3), nrow = 4, byrow = TRUE,
#'             dimnames = list(statesNames, statesNames))
#' mc <- new("markovchain", states = statesNames, transitionMatrix = P)
#' partition <- list(ab = c("a", "b"), cd = c("c", "d"))
#' is.lumpable(mc, partition) # rows of a and b give 0.3/0.7 to ab/cd
#' is.lumpable(mc, list(ac = c("a", "c"), bd = c("b", "d")))
#' @export
# Exact lumpability (Kemeny and Snell): the partition is lumpable if, for every group, the probability
# of moving into a group A is the same for all the states of the same group B.
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
#' @description Builds the chain on the blocks (macro-states) of a partition
#'   of the states. By default the chain must be exactly lumpable with respect
#'   to the partition (see [is.lumpable()]); with `force = TRUE` an
#'   approximate aggregation is returned for any partition.
#'
#' @param object A [`markovchain`][markovchain-class] object.
#' @param partition A named list of character vectors, the blocks of the
#'   partition; the names become the states of the result.
#' @param force Logical: if `FALSE` (default) stop unless the chain is
#'   exactly lumpable; if `TRUE` aggregate anyway.
#'
#' @details The probability of moving from block \eqn{A} to block \eqn{B} is
#'   the average over the states of \eqn{A} of their probability of moving
#'   into \eqn{B}, weighted by the stationary distribution (the average of the
#'   stationary distributions if there are several; equal weights within a
#'   block with zero stationary mass). For a lumpable chain all the states of
#'   a block give the same value, so the weights do not matter; otherwise the
#'   result is the usual stationary-weighted approximation, which preserves
#'   the stationary distribution of the blocks.
#'
#' @return A row-stochastic [`markovchain`][markovchain-class] object on the
#'   blocks.
#' @references Kemeny, J. G. and Snell, J. L. (1960). Finite Markov Chains.
#'   Van Nostrand.
#' @family state aggregation
#' @examples
#' statesNames <- c("a", "b", "c", "d")
#' P <- matrix(c(0.1, 0.2, 0.3, 0.4,
#'               0.2, 0.1, 0.4, 0.3,
#'               0.3, 0.3, 0.2, 0.2,
#'               0.2, 0.4, 0.1, 0.3), nrow = 4, byrow = TRUE,
#'             dimnames = list(statesNames, statesNames))
#' mc <- new("markovchain", states = statesNames, transitionMatrix = P)
#' partition <- list(ab = c("a", "b"), cd = c("c", "d"))
#' lump(mc, partition)
#' # an approximate aggregation over a partition that is not lumpable
#' lump(mc, list(ac = c("a", "c"), bd = c("b", "d")), force = TRUE)
#' @export
# lump(): builds the lumped chain. If the chain is not exactly lumpable, force = TRUE computes an
# approximate chain weighting the states by their stationary probability.
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

#' Aggregate a Markov chain by spectral clustering
#'
#' @description Finds a partition of the states into `k` blocks by k-means
#'   clustering of the states on the right eigenvectors of the `k`
#'   eigenvalues of largest modulus, then aggregates the chain over it with
#'   [`lump(force = TRUE)`][lump()]. States that the slow modes of the chain
#'   do not distinguish (metastable groups) end up in the same block.
#'
#' @param object A [`markovchain`][markovchain-class] object.
#' @param k Number of blocks, an integer between 2 and the number of states
#'   minus 1.
#'
#' @details This is a heuristic for approximate (metastable) aggregation, not
#'   a test of exact lumpability: check the result with [is.lumpable()], or
#'   use [aggregateStates()] for the aggregation that minimizes the
#'   Kullback-Leibler divergence rate. The k-means clustering is run with a
#'   fixed seed and 10 starts, so the result is reproducible; the random
#'   number stream of the session is left unchanged.
#'
#' @return A list with `partition`, a named list of the states of each block
#'   (`"Macro_1"`, `"Macro_2"`, ...), and `lumped_chain`, the aggregated
#'   [`markovchain`][markovchain-class] object.
#' @family state aggregation
#' @examples
#' # two groups of states that rarely communicate
#' statesNames <- c("a", "b", "c", "d")
#' P <- matrix(c(0.6, 0.38, 0.01, 0.01,
#'               0.3, 0.68, 0.01, 0.01,
#'               0.01, 0.01, 0.5, 0.48,
#'               0.01, 0.01, 0.4, 0.58), nrow = 4, byrow = TRUE,
#'             dimnames = list(statesNames, statesNames))
#' mc <- new("markovchain", states = statesNames, transitionMatrix = P)
#' autoLump(mc, k = 2)
#' @export
# autoLump(): finds a partition automatically by clustering (k-means) the states on the first k
# eigenvectors, then lumps the chain.
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

            # Eigenvectors of the k eigenvalues with the largest moduli.
            eig <- eigen(P)
            ord <- order(Mod(eig$values), decreasing = TRUE)
            V_eig <- Re(eig$vectors[, ord[seq_len(k)], drop = FALSE])

            # k-means starts at random: fix the seed (42) so that the result is reproducible, and restore the
            # user's random number stream when the function exits.
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

            # states with similar rows of the dominant eigenvectors go in the same block
            partition <- stats::setNames(vector("list", k), paste0("Macro_", seq_len(k)))
            for (i in seq_len(k)) {
              partition[[i]] <- state_names[clust$cluster == i]
            }

            # force = TRUE: the clustering gives an approximate lumping, not an exact one
            list(
              partition = partition,
              lumped_chain = lump(object, partition, force = TRUE)
            )
          })
