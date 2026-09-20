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
