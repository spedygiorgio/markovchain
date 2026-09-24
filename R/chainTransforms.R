#' Build a lazy version of a Markov chain
#'
#' Constructs the "lazy" chain associated with a \code{markovchain} object:
#' at every step, stay put with probability \code{alpha} and otherwise take
#' a step of the original chain.
#'
#' For a transition matrix \eqn{P} (in whichever storage convention
#' \code{object} already uses) and a laziness parameter
#' \eqn{\alpha\in[0,1]}, the lazy chain's transition matrix is
#' \deqn{L = \alpha I + (1-\alpha) P.}
#'
#' @param object A \code{markovchain} object.
#' @param alpha A single number in \eqn{[0,1]}, the probability of staying
#'   in the current state at each step. The default, \code{0.5}, matches
#'   the usual textbook "lazy random walk" construction.
#'
#' @return A new \code{markovchain} object with transition matrix
#'   \eqn{L = \alpha I + (1-\alpha) P}, the same states, and the same
#'   row/column-stochastic storage convention (\code{byrow}) as
#'   \code{object}.
#'
#' @details
#' Laziness is a standard device for forcing aperiodicity without changing
#' where the chain can go or its stationary distribution:
#' \itemize{
#'   \item \eqn{L} has the *same* stationary distribution as \eqn{P}
#'     (if \eqn{\pi P=\pi} then \eqn{\pi L = \alpha\pi + (1-\alpha)\pi P =
#'     \pi}), and the same communicating classes, since
#'     \eqn{L_{ij}>0 \iff P_{ij}>0} for \eqn{i\ne j}.
#'   \item For \eqn{0<\alpha<1}, \eqn{L} is aperiodic even if \eqn{P} is
#'     periodic, because \eqn{L_{ii}=\alpha>0} for every state \eqn{i}
#'     rules out any period greater than \eqn{1}. This is why
#'     \code{\link{mixingTime}}, which requires aperiodicity, is often
#'     applied to \code{lazyChain(object)} rather than to a periodic
#'     \code{object} directly (see \code{\link{mixingTime}}'s own
#'     documentation for why it rejects periodic chains outright rather
#'     than lazifying them automatically).
#'   \item Every non-trivial eigenvalue of \eqn{L} is
#'     \eqn{\alpha + (1-\alpha)\lambda} for the corresponding eigenvalue
#'     \eqn{\lambda} of \eqn{P}: laziness shrinks the whole non-trivial
#'     spectrum towards \eqn{\alpha}, so \code{\link{slem}} and
#'     \code{\link{impliedTimescales}} generally get *worse* (mixing gets
#'     slower) as \code{alpha} increases towards \code{1}.
#' }
#' \eqn{\alpha=0} returns \eqn{P} unchanged; \eqn{\alpha=1} returns the
#' identity matrix (a chain that never moves).
#'
#' @references
#' Levin, D. A. and Peres, Y. (2017). \emph{Markov Chains and Mixing Times},
#' 2nd edition. American Mathematical Society.
#'
#' @seealso \code{\link{subchain}}, \code{\link{mixingTime}},
#'   \code{\link{slem}}
#'
#' @examples
#' # A 2-cycle is periodic (period 2); its lazy version is aperiodic.
#' statesNames <- c("a", "b")
#' cycle2 <- new("markovchain", states = statesNames,
#'   transitionMatrix = matrix(c(0, 1, 1, 0), byrow = TRUE, nrow = 2,
#'                             dimnames = list(statesNames, statesNames)))
#' period(cycle2)
#' lazyCycle2 <- lazyChain(cycle2, alpha = 0.5)
#' period(lazyCycle2)
#' steadyStates(cycle2)
#' steadyStates(lazyCycle2) # unchanged by laziness
#'
#' @exportMethod lazyChain
setGeneric("lazyChain", function(object, alpha = 0.5) standardGeneric("lazyChain"))

#' @rdname lazyChain
setMethod("lazyChain", "markovchain", function(object, alpha = 0.5) {
  if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("alpha must be a single finite number in [0, 1].")
  }

  P <- as.matrix(object@transitionMatrix)
  n <- nrow(P)
  if (n != ncol(P) || any(!is.finite(P))) {
    stop("The transition matrix must be square and finite.")
  }

  # No byrow handling is needed here: I and P commute with transposition,
  # so alpha*I + (1-alpha)*P is row-stochastic exactly when P is, and
  # column-stochastic exactly when P is. The result inherits whichever
  # convention `object` already uses.
  L <- alpha * diag(n) + (1 - alpha) * P
  dimnames(L) <- dimnames(P)

  new("markovchain",
      states = states(object),
      byrow = object@byrow,
      transitionMatrix = L,
      name = paste0(object@name, " (lazy, alpha = ", alpha, ")"))
})

#' Restrict a Markov chain to a subset of states
#'
#' Restricts a \code{markovchain} object to a chosen subset of its states,
#' either as a raw (generally non-stochastic) principal submatrix, or as a
#' properly renormalized Markov chain describing behaviour conditional on
#' staying inside the subset.
#'
#' @param object A \code{markovchain} object.
#' @param states A character vector of at least one state name from
#'   \code{states(object)}, with no duplicates: the subset to restrict to.
#' @param method Either \code{"submatrix"} or \code{"renormalize"} (can be
#'   abbreviated). See Details: the two options answer genuinely different
#'   questions and are not interchangeable.
#'
#' @return
#' If \code{method = "submatrix"}: a plain numeric matrix (\emph{not} a
#' \code{markovchain} object, since its rows generally do not sum to one),
#' the principal submatrix of the transition matrix restricted to
#' \code{states}, always returned in row-stochastic orientation regardless
#' of \code{object}'s own storage convention.
#'
#' If \code{method = "renormalize"}: a new \code{markovchain} object on
#' exactly the states in \code{states}, row-stochastic, describing the
#' chain conditional on never leaving that subset.
#'
#' @details
#' The two methods are deliberately named after two different, standard
#' constructions, so that the choice -- and its consequences -- is explicit
#' rather than implied:
#'
#' \describe{
#'   \item{\code{"submatrix"}}{Simply the entries of \eqn{P} with both
#'     indices restricted to \code{states}, with no adjustment. Its rows
#'     generally sum to \emph{less than} one, because probability mass that
#'     originally went to states outside the subset is dropped, not
#'     redistributed. This is the "\eqn{Q}" block used, e.g., when building
#'     the fundamental matrix of an absorbing chain (see
#'     \code{\link{fundamentalMatrix}}): a useful building block for other
#'     computations, but not itself a transition matrix of any Markov
#'     chain, which is why it is returned as a plain matrix.}
#'   \item{\code{"renormalize"}}{Each retained row is divided by its own
#'     sum, so the result is row-stochastic and can be wrapped in a
#'     \code{markovchain} object. This is the chain of successive positions
#'     of \code{object}, \emph{conditioned on the event that it never
#'     leaves \code{states}} (sometimes called the chain "watched on"
#'     \code{states}, or its taboo probabilities; see Norris (1997),
#'     Section 3.3). It requires every state in \code{states} to have
#'     strictly positive probability of transitioning within the subset
#'     (otherwise that conditioning event has probability zero from that
#'     state, and the row cannot be renormalized); an error is raised
#'     naming any state that fails this, rather than silently producing a
#'     row of \code{NaN}.}
#' }
#'
#' Neither method requires \code{object} to be irreducible: restricting to
#' a subset of states is meaningful for any chain, and is often used
#' precisely to study one communicating class in isolation.
#'
#' The implementation performs no eigendecomposition; it is
#' \eqn{O(k^2)} time and memory for a subset of size \eqn{k}, after an
#' \eqn{O(n^2)} extraction from the full \eqn{n}-state matrix.
#'
#' @references
#' Norris, J. R. (1998). \emph{Markov Chains}. Cambridge University Press.
#'
#' @seealso \code{\link{lazyChain}}, \code{\link{fundamentalMatrix}},
#'   \code{\link{canonicForm}}
#'
#' @examples
#' statesNames <- c("a", "b", "c")
#' mc <- new("markovchain", states = statesNames,
#'   transitionMatrix = matrix(c(0.5, 0.3, 0.2,
#'                               0.2, 0.6, 0.2,
#'                               0.1, 0.1, 0.8), byrow = TRUE, nrow = 3,
#'                             dimnames = list(statesNames, statesNames)))
#'
#' # Raw submatrix: rows no longer sum to 1, mass has "leaked" to "c".
#' subchain(mc, c("a", "b"), method = "submatrix")
#' rowSums(subchain(mc, c("a", "b"), method = "submatrix"))
#'
#' # Renormalized: a genuine markovchain, conditional on staying in {a, b}.
#' watched <- subchain(mc, c("a", "b"), method = "renormalize")
#' watched
#' rowSums(watched@transitionMatrix)
#'
#' @exportMethod subchain
setGeneric("subchain", function(object, states, method = c("submatrix", "renormalize")) {
  standardGeneric("subchain")
})

#' @rdname subchain
setMethod("subchain", "markovchain", function(object, states, method = c("submatrix", "renormalize")) {
  method <- match.arg(method)

  if (!is.character(states) || length(states) < 1L || anyNA(states)) {
    stop("states must be a non-empty character vector with no missing values.")
  }
  if (anyDuplicated(states)) {
    stop("states must not contain duplicate state names.")
  }
  allStates <- states(object)
  missing <- setdiff(states, allStates)
  if (length(missing) > 0L) {
    stop("Unknown state(s): ", paste(missing, collapse = ", "))
  }

  P <- as.matrix(object@transitionMatrix)
  if (!object@byrow) {
    P <- t(P)
  }
  # From here on P is row-stochastic, regardless of object@byrow.

  sub <- P[states, states, drop = FALSE]
  dimnames(sub) <- list(states, states)

  if (method == "submatrix") {
    return(sub)
  }

  rowTotals <- rowSums(sub)
  bad <- states[rowTotals <= 0]
  if (length(bad) > 0L) {
    stop(paste0(
      "Cannot renormalize: state(s) ", paste(bad, collapse = ", "),
      " have zero probability of transitioning within the given subset, ",
      "so the chain conditioned on staying inside it is undefined there."
    ))
  }

  # Dividing the matrix by the vector `rowTotals` recycles it down each
  # column, i.e. row i is divided by rowTotals[i] -- the row-wise
  # renormalization we want, with no explicit diag() matrix needed.
  renormalized <- sub / rowTotals
  dimnames(renormalized) <- list(states, states)

  new("markovchain",
      states = states,
      byrow = TRUE,
      transitionMatrix = renormalized,
      name = paste0(object@name, " (subchain)"))
})

#' Merge two Markov chains by convex combination of their transition matrices
#'
#' Builds a new \code{markovchain} object whose transition matrix is a convex
#' combination of the transition matrices of two existing chains defined on
#' the same state space.
#'
#' For transition matrices \eqn{P_1} (from \code{object}) and \eqn{P_2}
#' (from \code{other}), and a blending factor \eqn{\gamma\in[0,1]}, the
#' merged transition matrix is
#' \deqn{P = (1-\gamma) P_1 + \gamma P_2.}
#'
#' @param object A \code{markovchain} object.
#' @param other A second \code{markovchain} object, defined on the same set
#'   of state names as \code{object} (see Details for what "same" means
#'   here).
#' @param gamma A single number in \eqn{[0,1]}: the weight given to
#'   \code{other}. \code{gamma = 0} returns \code{object} unchanged (up to
#'   storage convention); \code{gamma = 1} returns \code{other} unchanged.
#'
#' @return A new \code{markovchain} object, row-stochastic, on the common
#'   state names, with transition matrix \eqn{P = (1-\gamma)P_1+\gamma P_2}.
#'
#' @details
#' \strong{States are matched by name, not by position.} \code{object} and
#' \code{other} must have exactly the same set of state names (as sets --
#' \code{other}'s states may be in a different order, or \code{other} may use
#' a different storage convention (\code{byrow}), and both are handled
#' correctly). Rows and columns of \code{other}'s transition matrix are
#' realigned to \code{object}'s state order before combining, and both
#' matrices are converted to row-stochastic form first if needed, so that
#' \eqn{P_1} and \eqn{P_2} are always combined entry-for-entry between
#' matching states rather than between matching matrix positions.
#'
#' This is a deliberate difference from the na\"ive version of this
#' operation, which combines two same-\emph{size} transition matrices
#' positionally and would silently produce a meaningless result if the two
#' chains happened to list their states in a different order (or under a
#' different \code{byrow} convention) despite describing the same states.
#' Requiring identical state name sets, rather than merely identical size,
#' catches that mismatch as an error instead of propagating it.
#'
#' Because \eqn{P_1} and \eqn{P_2} are both row-stochastic with non-negative
#' entries and \eqn{\gamma\in[0,1]}, \eqn{P} is automatically row-stochastic
#' with non-negative entries: no renormalization is needed (this is the same
#' convexity argument used for \code{\link{lazyChain}}, of which
#' \code{mergeWith(object, identity_chain, gamma)} is a special case when
#' \code{other} is an identity chain on the same states).
#'
#' This function does not require \code{object} or \code{other} to be
#' irreducible: merging is meaningful for any two chains on the same state
#' space, including reducible ones. Note, however, that the merged chain's
#' stationary distribution (if any) is generally \emph{not} a combination of
#' \eqn{\pi_1} and \eqn{\pi_2} in any simple way; it must be recomputed from
#' \eqn{P} directly.
#'
#' The implementation performs no eigendecomposition and is \eqn{O(n^2)}
#' time and memory for two \eqn{n}-state chains, dominated by realigning
#' \code{other}'s matrix to \code{object}'s state order.
#'
#' @seealso \code{\link{lazyChain}}, \code{\link{subchain}}
#'
#' @examples
#' statesNames <- c("a", "b")
#' mc1 <- new("markovchain", states = statesNames,
#'   transitionMatrix = matrix(c(0.9, 0.1, 0.1, 0.9), byrow = TRUE, nrow = 2,
#'                             dimnames = list(statesNames, statesNames)))
#' # mc2 lists the same two states in the opposite order.
#' mc2 <- new("markovchain", states = rev(statesNames),
#'   transitionMatrix = matrix(c(0.5, 0.5, 0.5, 0.5), byrow = TRUE, nrow = 2,
#'                             dimnames = list(rev(statesNames), rev(statesNames))))
#' merged <- mergeWith(mc1, mc2, gamma = 0.5)
#' merged
#' # States are matched by name: merged["a", "b"] combines mc1["a","b"] with
#' # mc2["a","b"], not with whatever happened to sit in the same matrix cell.
#'
#' @exportMethod mergeWith
setGeneric("mergeWith", function(object, other, gamma = 0.5) standardGeneric("mergeWith"))

#' @rdname mergeWith
setMethod("mergeWith", signature(object = "markovchain", other = "markovchain"),
          function(object, other, gamma = 0.5) {
  if (length(gamma) != 1L || !is.numeric(gamma) || !is.finite(gamma) ||
      gamma < 0 || gamma > 1) {
    stop("gamma must be a single finite number in [0, 1].")
  }

  stateNames <- states(object)
  otherStateNames <- states(other)
  if (length(stateNames) != length(otherStateNames) ||
      !setequal(stateNames, otherStateNames)) {
    stop(paste0(
      "object and other must be defined on exactly the same set of state ",
      "names; found object states {", paste(sort(stateNames), collapse = ", "),
      "} and other states {", paste(sort(otherStateNames), collapse = ", "),
      "}."
    ))
  }

  P1 <- as.matrix(object@transitionMatrix)
  if (!object@byrow) {
    P1 <- t(P1)
  }
  P2 <- as.matrix(other@transitionMatrix)
  if (!other@byrow) {
    P2 <- t(P2)
  }
  # Realign other's matrix to object's state order by name, not position.
  P2 <- P2[stateNames, stateNames, drop = FALSE]

  n <- length(stateNames)
  if (any(!is.finite(P1)) || any(!is.finite(P2)) ||
      nrow(P1) != n || ncol(P1) != n) {
    stop("Both transition matrices must be square, finite, and match the number of states.")
  }

  P <- (1 - gamma) * P1 + gamma * P2
  dimnames(P) <- list(stateNames, stateNames)

  new("markovchain",
      states = stateNames,
      byrow = TRUE,
      transitionMatrix = P,
      name = paste0(object@name, " + ", other@name, " (merged, gamma = ", gamma, ")"))
})

#' Apply a boundary condition to a Markov chain's first and last state
#'
#' Replaces the transition rows of a \code{markovchain} object's first and
#' last state (in \code{states(object)} order) with an absorbing,
#' reflecting, or semi-reflecting rule, leaving every other row unchanged.
#'
#' @param object A \code{markovchain} object with at least 2 states.
#' @param boundaryCondition Either:
#'   \itemize{
#'     \item the string \code{"absorbing"}: the first and last state each
#'       become absorbing (\eqn{P_{11}=1}, \eqn{P_{nn}=1});
#'     \item the string \code{"reflecting"}: the first state moves to the
#'       second with certainty and the last state moves to the
#'       second-to-last with certainty (\eqn{P_{12}=1},
#'       \eqn{P_{n,n-1}=1});
#'     \item a single number \eqn{\beta\in[0,1]}, the
#'       \emph{semi-reflecting} case: the first state stays with
#'       probability \eqn{1-\beta} and moves to the second state with
#'       probability \eqn{\beta} (\eqn{P_{11}=1-\beta}, \eqn{P_{12}=\beta}),
#'       and symmetrically the last state stays with probability
#'       \eqn{1-\beta} and moves to the second-to-last with probability
#'       \eqn{\beta}. \eqn{\beta=0} is the absorbing case and \eqn{\beta=1}
#'       is the reflecting case.
#'   }
#'
#' @return A new \code{markovchain} object, row-stochastic, on the same
#'   states as \code{object}, identical to \code{object} except in its
#'   first and last transition rows.
#'
#' @details
#' This function assumes -- as is standard for a boundary condition -- that
#' \code{states(object)} is meaningfully ordered along a line, first state
#' to last state, as it would be e.g. for \code{\link{birthDeath}} or any
#' other chain built to represent a bounded random walk. It does not check
#' this (there is no general way to check it from the transition matrix
#' alone) and applies the same first/last-row replacement regardless of
#' \code{object}'s actual structure; only the two boundary rows are ever
#' touched, so applying it to a chain whose states are not linearly ordered
#' simply reinterprets whichever states happen to be listed first and last.
#'
#' Unlike \code{\link{gamblersRuin}}, which is absorbing at both ends by
#' construction and cannot be un-done, \code{toBoundedChain()} can be
#' applied to any existing chain and with any of the three conditions,
#' including reflecting or semi-reflecting ones that \code{gamblersRuin()}
#' does not offer directly.
#'
#' The implementation touches only 2 of the \eqn{n} rows and is
#' \eqn{O(n)} time and memory beyond copying the transition matrix.
#'
#' @seealso \code{\link{birthDeath}}, \code{\link{gamblersRuin}}
#'
#' @examples
#' bd <- birthDeath(p = c(0.3, 0.4, 0.5), q = c(0.2, 0.3, 0.1))
#'
#' absorbed <- toBoundedChain(bd, "absorbing")
#' absorbed@transitionMatrix[1, ]
#' absorbed@transitionMatrix[4, ]
#'
#' reflected <- toBoundedChain(bd, "reflecting")
#' reflected@transitionMatrix[1, ]
#'
#' semiReflected <- toBoundedChain(bd, 0.25)
#' semiReflected@transitionMatrix[1, ]
#'
#' @exportMethod toBoundedChain
setGeneric("toBoundedChain", function(object, boundaryCondition) {
  standardGeneric("toBoundedChain")
})

#' @rdname toBoundedChain
setMethod("toBoundedChain", "markovchain", function(object, boundaryCondition) {
  stateNames <- states(object)
  n <- length(stateNames)
  if (n < 2L) {
    stop("toBoundedChain requires a chain with at least 2 states.")
  }

  P <- as.matrix(object@transitionMatrix)
  if (!object@byrow) {
    P <- t(P)
  }

  firstRow <- numeric(n)
  lastRow <- numeric(n)

  if (is.character(boundaryCondition) && length(boundaryCondition) == 1L &&
      boundaryCondition %in% c("absorbing", "reflecting")) {
    if (boundaryCondition == "absorbing") {
      firstRow[1] <- 1
      lastRow[n] <- 1
    } else {
      firstRow[min(2L, n)] <- 1
      lastRow[max(n - 1L, 1L)] <- 1
    }
  } else if (is.numeric(boundaryCondition) && length(boundaryCondition) == 1L &&
             is.finite(boundaryCondition) &&
             boundaryCondition >= 0 && boundaryCondition <= 1) {
    beta <- boundaryCondition
    firstRow[1] <- 1 - beta
    firstRow[min(2L, n)] <- firstRow[min(2L, n)] + beta
    lastRow[n] <- 1 - beta
    lastRow[max(n - 1L, 1L)] <- lastRow[max(n - 1L, 1L)] + beta
  } else {
    stop(paste0(
      "boundaryCondition must be \"absorbing\", \"reflecting\", or a ",
      "single number in [0, 1]."
    ))
  }

  P[1, ] <- firstRow
  P[n, ] <- lastRow
  dimnames(P) <- list(stateNames, stateNames)

  new("markovchain",
      states = stateNames,
      byrow = TRUE,
      transitionMatrix = P,
      name = paste0(object@name, " (bounded)"))
})

#' Return the n-step transition chain
#'
#' Returns the \code{markovchain} object whose transition matrix is
#' \eqn{P^{\code{order}}}: from any state, its one-step transition
#' probabilities are the original chain's \code{order}-step transition
#' probabilities.
#'
#' @param object A \code{markovchain} object.
#' @param order A single integer of at least \code{2}.
#'
#' @return A new \code{markovchain} object on the same states as
#'   \code{object}, with transition matrix \eqn{P^{\code{order}}}.
#'
#' @details
#' This is a thin, discoverability-only wrapper around
#' \code{object ^ order} (see \code{\link[=markovchain-class]{^,markovchain,numeric-method}}),
#' provided under this name because \pkg{PyDTMC}'s equivalent method is
#' called \code{to_nth_order()}. It exists so that the operation is easy to
#' find by that name; it introduces no new computation; the underlying
#' \code{^} method is already \eqn{O(n^3\log(\code{order}))} via repeated
#' squaring (\code{expm::\link[expm]{\%^\%}}), not a naive \code{order}-fold
#' product, so there is nothing to improve on algorithmically here.
#'
#' @seealso \code{\link{toBoundedChain}}, \code{\link{lazyChain}}
#'
#' @examples
#' mc <- new("markovchain", states = c("a", "b"),
#'           transitionMatrix = matrix(c(0.9, 0.1, 0.3, 0.7), byrow = TRUE, nrow = 2))
#' identical(unclass(toNthOrder(mc, 5)@transitionMatrix), unclass((mc ^ 5)@transitionMatrix))
#'
#' @export
toNthOrder <- function(object, order) {
  if (!is(object, "markovchain")) {
    stop("object must be a markovchain object.")
  }
  if (length(order) != 1L || !is.numeric(order) || is.na(order) ||
      order != as.integer(order) || order < 2L) {
    stop("order must be a single integer of at least 2.")
  }
  object ^ as.integer(order)
}
