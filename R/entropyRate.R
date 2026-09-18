#' Entropy rate of a Markov chain
#'
#' Computes the entropy rate of a finite, irreducible discrete-time Markov
#' chain from its stationary distribution and transition matrix.
#'
#' For a row-stochastic transition matrix \eqn{P} and stationary distribution
#' \eqn{\pi}, the entropy rate is
#' \deqn{H = -\sum_i \pi_i \sum_j p_{ij}\log_b(p_{ij}),}
#' with zero-probability transitions contributing zero by continuity.
#'
#' @param object A \code{markovchain} object representing a finite,
#'   irreducible discrete-time Markov chain.
#' @param base A finite numeric scalar strictly greater than one. The default,
#'   \code{2}, returns entropy in bits per transition. Use \code{exp(1)} for
#'   nats per transition.
#'
#' @return A non-negative numeric scalar containing the entropy rate in units
#'   determined by \code{base}.
#'
#' @details
#' For a stationary first-order Markov chain, the entropy rate equals the
#' conditional entropy \eqn{H(X_{t+1}\mid X_t)}. Irreducibility guarantees a
#' unique stationary distribution; aperiodicity is not required.
#'
#' Reducible chains may admit multiple stationary distributions and therefore
#' different entropy rates. This method rejects them rather than silently
#' selecting one stationary distribution.
#'
#' Transitions with probability zero are ignored, implementing the standard
#' convention \eqn{0\log(0)=0} without evaluating \code{log(0)}.
#'
#' The stationary-distribution computation dominates the running time. Once
#' the stationary distribution is available, evaluating the entropy rate takes
#' \eqn{O(n^2)} time and \eqn{O(n^2)} temporary memory for a dense \eqn{n}-state
#' transition matrix.
#'
#' @references
#' Cover, T. M. and Thomas, J. A. (2006). \emph{Elements of Information
#' Theory}, 2nd edition. Wiley.
#'
#' Strelioff, C. C., Crutchfield, J. P. and Huebler, A. W. (2007).
#' Inferring Markov chains: Bayesian estimation, model comparison, entropy
#' rate, and out-of-class modeling. \emph{Physical Review E}, 76, 011106.
#'
#' @seealso \code{\link{steadyStates}}, \code{\link{is.irreducible}}
#'
#' @examples
#' statesNames <- c("a", "b")
#' mc <- new("markovchain",
#'   states = statesNames,
#'   transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
#'     byrow = TRUE, nrow = 2,
#'     dimnames = list(statesNames, statesNames)))
#' entropyRate(mc)
#' entropyRate(mc, base = exp(1))
#'
#' @exportMethod entropyRate
setGeneric("entropyRate", function(object, base = 2) standardGeneric("entropyRate"))

#' @rdname entropyRate
setMethod("entropyRate", "markovchain", function(object, base = 2) {
  if (length(base) != 1L || !is.numeric(base) || !is.finite(base) ||
      base <= 1) {
    stop("base must be a single finite number strictly greater than 1.")
  }
  if (!is.irreducible(object)) {
    stop(paste0(
      "Entropy rate is defined here using the unique stationary ",
      "distribution of an irreducible Markov chain."
    ))
  }

  # All computations below use a row-stochastic representation. The state
  # order is unchanged by this transpose, so it remains aligned with pi.
  P <- as.matrix(object@transitionMatrix)
  if (!object@byrow) {
    P <- t(P)
  }
  n <- nrow(P)
  if (n != ncol(P) || any(!is.finite(P)) || any(P < 0)) {
    stop("The transition matrix must be square, finite, and non-negative.")
  }

  pi <- as.numeric(steadyStates(object))
  if (length(pi) != n || any(!is.finite(pi)) ||
      any(pi < -sqrt(.Machine$double.eps)) || sum(pi) <= 0) {
    stop("Unable to obtain a valid stationary distribution.")
  }
  # Remove harmless negative round-off and normalize defensively.
  pi[pi < 0] <- 0
  pi <- pi / sum(pi)

  # Store p_ij log_b(p_ij) directly, avoiding separate full-size log(P) and
  # P * log(P) matrices. Entries with p_ij = 0 remain zero.
  contributions <- matrix(0, nrow = n, ncol = n)
  positive <- P > 0
  logBase <- log(base)
  contributions[positive] <-
    P[positive] * log(P[positive]) / logBase

  value <- -sum(pi * rowSums(contributions))
  # Guard only against a negative signed zero or round-off at zero.
  if (value < 0 && value > -100 * .Machine$double.eps) {
    value <- 0
  }
  if (!is.finite(value) || value < 0) {
    stop("Unable to compute a finite non-negative entropy rate.")
  }
  as.numeric(value)
})
