#' Fundamental matrix of an absorbing Markov chain
#'
#' Computes the fundamental matrix of a finite absorbing discrete-time Markov
#' chain. If \eqn{Q} is the transition submatrix restricted to transient
#' states, the fundamental matrix is
#' \deqn{N = I + Q + Q^2 + \cdots = (I - Q)^{-1}.}
#' The \eqn{(i,j)} entry is the expected number of visits to transient state
#' \eqn{j}, including the initial visit when \eqn{i = j}, before absorption,
#' when the chain starts in transient state \eqn{i}.
#'
#' @param object A \code{markovchain} object representing an absorbing
#'   Markov chain.
#'
#' @return A numeric matrix containing the fundamental matrix, with transient
#'   state names as row and column names. If all states are absorbing, the
#'   result is a 0-by-0 matrix because there are no transient states.
#'
#' @details
#' For a finite absorbing chain, the state space can be reordered so that the
#' transition matrix has canonical form with transient block \eqn{Q} and an
#' absorbing block. The spectral radius of \eqn{Q} is less than one, so the
#' Neumann series converges and \eqn{I-Q} is nonsingular.
#'
#' The fundamental matrix also gives the expected time to absorption through
#' \eqn{t = N 1}, where \eqn{1} is a vector of ones, and absorption
#' probabilities through \eqn{B = N R}, where \eqn{R} contains transition
#' probabilities from transient to absorbing states.
#'
#' The function requires an absorbing Markov chain: at least one absorbing
#' state must exist and every recurrent state must be absorbing. A chain with
#' no transient states is a valid degenerate case and returns a 0-by-0 matrix.
#'
#' @references
#' Kemeny, J. G. and Snell, J. L. (1976). *Finite Markov Chains*.
#' Springer.
#'
#' Grinstead, C. M. and Snell, J. L. (1997). *Introduction to Probability*.
#' American Mathematical Society.
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
