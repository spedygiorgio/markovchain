#' Fundamental matrix of an absorbing Markov chain
#'
#' Computes the fundamental matrix of an absorbing discrete-time Markov chain.
#' For the transient-state submatrix \eqn{Q}, the fundamental matrix is
#' \eqn{N = (I - Q)^{-1}}. Its \eqn{(i,j)} entry is the expected number of
#' visits to transient state \eqn{j} before absorption when starting from
#' transient state \eqn{i}.
#'
#' @param object A \code{markovchain} object representing an absorbing
#'   Markov chain.
#'
#' @return A numeric matrix containing the fundamental matrix, with transient
#'   state names as row and column names.
#'
#' @details
#' An absorbing Markov chain must contain at least one absorbing state, and all
#' recurrent states must be absorbing. The function extracts the transition
#' submatrix corresponding to transient states and computes its inverse
#' complement. If there are no transient states, an empty matrix is returned.
#'
#' @references
#' Kemeny, J. G. and Snell, J. L. (1976). *Finite Markov Chains*.
#' Springer.
#'
#' @seealso \code{\link{absorbingStates}}, \code{\link{transientStates}},
#'   \code{\link{meanAbsorptionTime}}, \code{\link{absorptionProbabilities}}
#'
#' @examples
#' states <- c("a", "b", "absorbed")
#' mc <- new("markovchain", states = states,
#'   transitionMatrix = matrix(c(
#'     0.5, 0.4, 0.1,
#'     0.2, 0.6, 0.2,
#'     0,   0,   1
#'   ), nrow = 3, byrow = TRUE,
#'   dimnames = list(states, states)))
#'
#' fundamentalMatrix(mc)
#'
#' @export
fundamentalMatrix <- function(object) {
  if (!is(object, "markovchain")) {
    stop("please provide a valid markovchain object")
  }

  absorbing <- absorbingStates(object)
  recurrent <- recurrentStates(object)

  if (length(absorbing) == 0L || !all(recurrent %in% absorbing)) {
    stop("fundamental matrix requires an absorbing Markov chain")
  }

  transient <- transientStates(object)
  if (length(transient) == 0L) {
    out <- matrix(numeric(0), nrow = 0L, ncol = 0L)
    dimnames(out) <- list(character(0), character(0))
    return(out)
  }

  q <- object@transitionMatrix[transient, transient, drop = FALSE]
  out <- solve(diag(nrow(q)) - q)
  dimnames(out) <- list(transient, transient)
  out
}
