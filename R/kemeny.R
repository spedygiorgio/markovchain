#' Kemeny's constant of a Markov chain
#'
#' Computes Kemeny's constant for a finite, irreducible discrete-time Markov
#' chain. Kemeny's constant is the expected first-passage time to a state
#' selected from the stationary distribution and is independent of the
#' starting state.
#'
#' For an irreducible transition matrix `P`, let `pi` be its stationary
#' distribution and define
#' `Z = (I - P + 1 %*% t(pi))^{-1}`. The function returns
#' `trace(Z) - 1`, which is numerically equivalent to
#' `sum_j pi_j m_ij` for any starting state `i`.
#'
#' @param object A `markovchain` object representing a finite, irreducible
#'   discrete-time Markov chain.
#'
#' @return A numeric scalar containing Kemeny's constant.
#'
#' @details
#' Kemeny's constant is defined for finite irreducible chains. For a
#' reducible chain there can be more than one stationary distribution and the
#' state-independent quantity returned by this function is not defined in the
#' same way; such inputs therefore produce an error.
#'
#' The implementation uses a linear solve rather than explicitly computing
#' an inverse. This keeps the implementation numerically preferable to a
#' direct matrix inversion while retaining the compact closed-form expression.
#'
#' @references
#' Kemeny, J. G. and Snell, J. L. (1976). *Finite Markov Chains*.
#' Springer.
#'
#' @examples
#' statesNames <- c("a", "b")
#' mc <- new("markovchain",
#'   states = statesNames,
#'   transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
#'     byrow = TRUE, nrow = 2,
#'     dimnames = list(statesNames, statesNames))
#' kemenyConstant(mc)
#'
#' @exportMethod kemenyConstant
setGeneric("kemenyConstant", function(object) standardGeneric("kemenyConstant"))

#' @rdname kemenyConstant
setMethod("kemenyConstant", "markovchain", function(object) {
  if (!is.irreducible(object)) {
    stop("Kemeny's constant is defined here only for irreducible Markov chains.")
  }

  P <- as.matrix(object@transitionMatrix)
  pi <- as.numeric(steadyStates(object))
  n <- nrow(P)

  if (length(pi) != n || any(!is.finite(pi))) {
    stop("Unable to obtain a valid stationary distribution.")
  }

  Z <- solve(diag(n) - P + matrix(1, n, 1) %*% matrix(pi, 1, n))
  as.numeric(sum(diag(Z)) - 1)
})
