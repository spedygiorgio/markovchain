#' Entropy rate of a Markov chain
#'
#' Computes the entropy rate of a finite, irreducible discrete-time Markov
#' chain from its stationary distribution and transition matrix.
#'
#' For stationary distribution `pi` and transition matrix `P`, the entropy
#' rate is
#' `H = -sum_i pi_i sum_j p_ij log(p_ij)`, with zero-probability transitions
#' contributing zero by continuity.
#'
#' @param object A `markovchain` object representing a finite, irreducible
#'   discrete-time Markov chain.
#' @param base A positive numeric scalar specifying the logarithm base.
#'   The default, `2`, returns entropy in bits. Use `exp(1)` for nats.
#'
#' @return A numeric scalar containing the entropy rate in units determined by
#'   `base`.
#'
#' @details
#' The implementation uses the unique stationary distribution of an
#' irreducible finite chain. For reducible chains, the entropy rate depends on
#' which stationary distribution is selected, so such inputs produce an error
#' rather than silently choosing one of several stationary distributions.
#'
#' Transitions with probability zero are ignored, avoiding `log(0)` warnings
#' and implementing the standard convention `0 log(0) = 0`.
#'
#' @references
#' Cover, T. M. and Thomas, J. A. (2006). *Elements of Information Theory*.
#' Wiley.
#'
#' @examples
#' statesNames <- c("a", "b")
#' mc <- new("markovchain",
#'   states = statesNames,
#'   transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
#'     byrow = TRUE, nrow = 2,
#'     dimnames = list(statesNames, statesNames))
#' entropyRate(mc)
#' entropyRate(mc, base = exp(1))
#'
#' @exportMethod entropyRate
setGeneric("entropyRate", function(object, base = 2) standardGeneric("entropyRate"))

#' @rdname entropyRate
setMethod("entropyRate", "markovchain", function(object, base = 2) {
  if (length(base) != 1L || !is.numeric(base) || !is.finite(base) ||
      base <= 0 || base == 1) {
    stop("base must be a single finite positive number different from 1.")
  }
  if (!is.irreducible(object)) {
    stop("Entropy rate is defined here using the unique stationary distribution of irreducible Markov chains.")
  }

  P <- as.matrix(object@transitionMatrix)
  pi <- as.numeric(steadyStates(object))

  if (length(pi) != nrow(P) || any(!is.finite(pi))) {
    stop("Unable to obtain a valid stationary distribution.")
  }

  logP <- matrix(0, nrow(P), ncol(P))
  positive <- P > 0
  logP[positive] <- log(P[positive], base = base)
  stateEntropy <- -rowSums(P * logP)
  as.numeric(sum(pi * stateEntropy))
})
