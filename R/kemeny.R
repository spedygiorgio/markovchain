#' Kemeny's constant of a Markov chain
#'
#' Computes Kemeny's constant for a finite, irreducible discrete-time Markov
#' chain. It is the stationary-distribution-weighted mean hitting time of a
#' randomly selected destination and is independent of the starting state.
#'
#' For a row-stochastic transition matrix \eqn{P}, let \eqn{\pi} be its unique
#' stationary distribution and define
#' \deqn{Z = (I - P + \mathbf{1}\pi^T)^{-1}.}
#' With hitting times defined by
#' \eqn{T_j = \inf\{n \ge 0: X_n=j\}}, so that \eqn{m_{jj}=0}, the function
#' returns
#' \deqn{K = \sum_j \pi_j m_{ij} = \mathrm{tr}(Z)-1.}
#' The value does not depend on the starting state \eqn{i}.
#'
#' @param object A \code{markovchain} object representing a finite, irreducible
#'   discrete-time Markov chain.
#'
#' @return A numeric scalar containing Kemeny's constant.
#'
#' @details
#' Irreducibility is sufficient; aperiodicity is not required. Reducible chains
#' can have multiple stationary distributions and are rejected.
#'
#' Some references instead put the mean first-return time
#' \eqn{m_{jj}=1/\pi_j} on the diagonal. Under that convention the corresponding
#' stationary weighted sum is \eqn{K+1}, not \eqn{K}. This function uses the
#' zero-diagonal hitting-time convention, consistently with
#' \code{meanFirstPassageTime()}.
#'
#' The implementation uses a dense LAPACK solve for the fundamental matrix
#' \eqn{Z}. Its time complexity is \eqn{O(n^3)} and its memory use is
#' \eqn{O(n^2)}, as expected for a dense exact computation. It supports both
#' row- and column-stochastic storage.
#'
#' @references
#' Kemeny, J. G. and Snell, J. L. (1960). \emph{Finite Markov Chains}.
#' D. Van Nostrand, Princeton, NJ.
#'
#' @seealso \code{\link{meanFirstPassageTime}},
#'   \code{\link{steadyStates}}, \code{\link{is.irreducible}}
#'
#' @examples
#' statesNames <- c("a", "b")
#' mc <- new("markovchain",
#'   states = statesNames,
#'   transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
#'     byrow = TRUE, nrow = 2,
#'     dimnames = list(statesNames, statesNames)))
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
  if (!object@byrow) {
    P <- t(P)
  }
  n <- nrow(P)
  if (n != ncol(P) || any(!is.finite(P))) {
    stop("The transition matrix must be square and finite.")
  }

  pi <- as.numeric(steadyStates(object))
  if (length(pi) != n || any(!is.finite(pi)) || sum(pi) <= 0) {
    stop("Unable to obtain a valid stationary distribution.")
  }
  pi <- pi / sum(pi)

  # Add pi_j to every entry of column j without explicitly constructing
  # the n-by-n outer product 1 %*% t(pi).
  A <- diag(n) - P
  A <- sweep(A, 2L, pi, FUN = "+")
  Z <- solve(A)
  value <- as.numeric(sum(diag(Z)) - 1)
  if (!is.finite(value)) {
    stop("Unable to compute a finite Kemeny constant.")
  }
  value
})
