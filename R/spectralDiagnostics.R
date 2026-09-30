# Internal helper shared by slem(), spectralGap() and impliedTimescales().
#
# For a finite, irreducible, row-stochastic transition matrix P, the
# eigenvalue 1 has algebraic multiplicity exactly one (Perron-Frobenius).
# This holds regardless of periodicity: a periodic irreducible chain still
# has a *single* eigenvalue equal to 1, even though it additionally has
# other eigenvalues on the unit circle (e.g. -1 for a 2-cycle). This helper
# therefore removes only the one eigenvalue closest to 1 -- not every
# eigenvalue of unit modulus -- so that periodic chains are handled
# correctly instead of being rejected outright.
#
# Returns the moduli of the remaining n-1 eigenvalues, sorted in decreasing
# order. For a one-state chain it returns numeric(0).
.nonTrivialEigenvalueModuli <- function(object, caller) {
  if (!is.irreducible(object)) {
    stop(paste0(
      caller, " is defined here only for irreducible Markov chains."
    ))
  }

  P <- as.matrix(object@transitionMatrix)
  if (!object@byrow) {
    P <- t(P)
  }
  n <- nrow(P)
  if (n != ncol(P) || any(!is.finite(P))) {
    stop("The transition matrix must be square and finite.")
  }

  if (n == 1L) {
    return(numeric(0))
  }

  # Eigenvectors are not needed by any of the callers.
  values <- eigen(P, only.values = TRUE)$values

  unit <- which.min(Mod(values - 1))
  tol <- sqrt(.Machine$double.eps) * max(1, n)
  if (Mod(values[unit] - 1) > tol) {
    stop("Unable to identify the trivial (unit) eigenvalue of the transition matrix.")
  }

  moduli <- Mod(values[-unit])
  sort(moduli, decreasing = TRUE)
}

#' Second largest eigenvalue modulus (SLEM) of a Markov chain
#'
#' Computes the second largest eigenvalue modulus (SLEM) of a finite,
#' irreducible discrete-time Markov chain.
#'
#' For a row-stochastic transition matrix \eqn{P}, let
#' \eqn{1=\lambda_1,\lambda_2,\ldots,\lambda_n} be its eigenvalues. By the
#' Perron-Frobenius theorem an irreducible chain has \eqn{\lambda_1=1} with
#' algebraic multiplicity one, and \eqn{|\lambda_k|\le 1} for every
#' \eqn{k}. The SLEM is
#' \deqn{\mathrm{SLEM} = \max_{k>1} |\lambda_k|.}
#'
#' @param object A \code{markovchain} object representing a finite,
#'   irreducible discrete-time Markov chain.
#'
#' @return A numeric scalar in \eqn{[0,1]} containing the SLEM. For the
#'   trivial one-state chain, \code{0} is returned.
#'
#' @details
#' Only irreducibility is required, not aperiodicity. If the chain is
#' periodic, at least one non-trivial eigenvalue also has modulus one (e.g.
#' \eqn{\lambda=-1} for a 2-cycle), so \code{slem()} correctly returns
#' \code{1} rather than rejecting the chain: a periodic chain genuinely does
#' not contract towards its stationary distribution, which \code{SLEM = 1}
#' reflects.
#'
#' Repeated or complex non-trivial eigenvalues are handled through their
#' modulus \code{Mod()}, so complex-conjugate pairs contribute the same
#' value and ties do not need to be broken.
#'
#' The implementation calls \code{eigen()} with \code{only.values = TRUE},
#' so it never computes eigenvectors. Its time complexity is
#' \eqn{O(n^3)} and its memory use is \eqn{O(n^2)} for a dense \eqn{n}-state
#' transition matrix. It supports both row- and column-stochastic storage.
#'
#' @references
#' Levin, D. A. and Peres, Y. (2017). \emph{Markov Chains and Mixing Times},
#' 2nd edition. American Mathematical Society.
#'
#' @seealso \code{\link{spectralGap}}, \code{\link{impliedTimescales}},
#'   \code{\link{is.irreducible}}, \code{\link{period}}
#'
#' @examples
#' statesNames <- c("a", "b")
#' mc <- new("markovchain",
#'   states = statesNames,
#'   transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
#'     byrow = TRUE, nrow = 2,
#'     dimnames = list(statesNames, statesNames)))
#' slem(mc)
#'
#' @exportMethod slem
setGeneric("slem", function(object) standardGeneric("slem"))

#' @rdname slem
setMethod("slem", "markovchain", function(object) {
  moduli <- .nonTrivialEigenvalueModuli(object, "slem")
  value <- if (length(moduli) == 0L) 0 else moduli[1]
  if (!is.finite(value)) {
    stop("Unable to compute a finite SLEM.")
  }
  as.numeric(value)
})

#' Spectral gap of a Markov chain
#'
#' Computes the spectral gap of a finite, irreducible discrete-time Markov
#' chain, a lightweight diagnostic of its convergence and mixing behaviour.
#'
#' The spectral gap is defined from the second largest eigenvalue modulus
#' (SLEM, see \code{\link{slem}}) as
#' \deqn{\mathrm{gap} = 1 - \mathrm{SLEM}.}
#'
#' @param object A \code{markovchain} object representing a finite,
#'   irreducible discrete-time Markov chain.
#'
#' @return A numeric scalar in \eqn{[0,1]} containing the spectral gap. For
#'   the trivial one-state chain, \code{1} is returned.
#'
#' @details
#' As with \code{\link{slem}}, only irreducibility is required. A periodic
#' chain has \code{SLEM = 1} and therefore spectral gap \code{0}: this is
#' the mathematically correct value, not an error condition, since a
#' periodic chain never contracts towards its stationary distribution.
#' A larger spectral gap indicates faster convergence to stationarity; see
#' \code{\link{impliedTimescales}} for the timescale associated with each
#' non-trivial eigenvalue individually, of which the SLEM gives the slowest
#' (dominant) one.
#'
#' @references
#' Levin, D. A. and Peres, Y. (2017). \emph{Markov Chains and Mixing Times},
#' 2nd edition. American Mathematical Society.
#'
#' @seealso \code{\link{slem}}, \code{\link{impliedTimescales}},
#'   \code{\link{is.irreducible}}
#'
#' @examples
#' statesNames <- c("a", "b")
#' mc <- new("markovchain",
#'   states = statesNames,
#'   transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
#'     byrow = TRUE, nrow = 2,
#'     dimnames = list(statesNames, statesNames)))
#' spectralGap(mc)
#'
#' @exportMethod spectralGap
setGeneric("spectralGap", function(object) standardGeneric("spectralGap"))

#' @rdname spectralGap
setMethod("spectralGap", "markovchain", function(object) {
  as.numeric(1 - slem(object))
})

#' Implied timescales of a Markov chain
#'
#' Computes the implied relaxation timescale associated with each
#' non-trivial eigenvalue of a finite, irreducible discrete-time Markov
#' chain.
#'
#' For a row-stochastic transition matrix \eqn{P} with eigenvalues
#' \eqn{1=\lambda_1,\lambda_2,\ldots,\lambda_n} (\eqn{|\lambda_1|} the unique
#' unit eigenvalue of an irreducible chain), the implied timescale of
#' \eqn{\lambda_k}, \eqn{k>1}, is
#' \deqn{\tau_k = -\frac{1}{\log|\lambda_k|}}
#' for \eqn{0<|\lambda_k|<1}. Each \eqn{\tau_k} measures how many steps the
#' mode associated with \eqn{\lambda_k} takes to decay by a factor of
#' \eqn{1/e}; larger timescales correspond to slower-decaying, more
#' persistent modes.
#'
#' @param object A \code{markovchain} object representing a finite,
#'   irreducible discrete-time Markov chain.
#'
#' @return A named numeric vector of length \eqn{n-1} (one entry per
#'   non-trivial eigenvalue), sorted by decreasing timescale, i.e. by
#'   decreasing eigenvalue modulus. Names are \code{"tau2"}, \code{"tau3"},
#'   \ldots, matching the usual eigenvalue indexing
#'   \eqn{\lambda_2,\lambda_3,\ldots} in decreasing modulus. For the trivial
#'   one-state chain, a length-zero named numeric vector is returned.
#'
#' @details
#' Only irreducibility is required, not aperiodicity: this is the same
#' convention used by \code{\link{slem}} and \code{\link{spectralGap}}, and
#' it lets \code{impliedTimescales()} document periodic and boundary cases
#' explicitly rather than rejecting them:
#' \itemize{
#'   \item If \eqn{|\lambda_k|} is (numerically) exactly \eqn{1} -- which
#'     happens for non-trivial eigenvalues of periodic chains, e.g.
#'     \eqn{\lambda=-1} for a 2-cycle -- the corresponding mode never
#'     decays and \code{tau_k = Inf} is returned. This is a boundary case
#'     of the formula above (as \eqn{|\lambda|\to 1^-}, \eqn{\tau\to\infty})
#'     that is handled explicitly rather than by evaluating
#'     \eqn{-1/\log(1)}, which is numerically \code{-Inf} rather than the
#'     mathematically correct \code{+Inf}.
#'   \item If \eqn{|\lambda_k|} is (numerically) exactly \eqn{0}, the mode
#'     decays immediately and \code{tau_k = 0} is returned. \code{log(0)}
#'     evaluates to \code{-Inf} in \R, so this case is already handled
#'     correctly by the formula itself and needs no special-casing.
#' }
#'
#' The term "implied timescale" follows the Markov state model literature
#' in molecular kinetics, where it is additionally used, across chains
#' estimated at increasing lag times, as a self-consistency check on the
#' Markov (memoryless) approximation: implied timescales that are
#' approximately constant across lag times support the model, while ones
#' that drift indicate it should be revisited (see Prinz et al. (2011)).
#' Building such a lag-time comparison is left to the user, since it
#' requires re-estimating the chain at each lag: \code{impliedTimescales()}
#' itself only evaluates a single, already-fitted \code{markovchain} object.
#'
#' The implementation calls \code{eigen()} with \code{only.values = TRUE},
#' so it never computes eigenvectors. Its time complexity is
#' \eqn{O(n^3)} and its memory use is \eqn{O(n^2)} for a dense \eqn{n}-state
#' transition matrix. It supports both row- and column-stochastic storage.
#'
#' @references
#' Swope, W. C., Pitera, J. W. and Suits, F. (2004). Describing protein
#' folding kinetics by molecular dynamics simulations, 1: Theory.
#' \emph{J. Phys. Chem. B}, 108(21), 6571-6581.
#'
#' Prinz, J.-H., Wu, H., Sarich, M., Keller, B., Senne, M., Held, M.,
#' Chodera, J. D., Schutte, C. and Noe, F. (2011). Markov models of
#' molecular kinetics: Generation and validation. \emph{Journal of Chemical
#' Physics}, 134(17), 174105.
#'
#' @seealso \code{\link{slem}}, \code{\link{spectralGap}},
#'   \code{\link{is.irreducible}}
#'
#' @examples
#' statesNames <- c("a", "b")
#' mc <- new("markovchain",
#'   states = statesNames,
#'   transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
#'     byrow = TRUE, nrow = 2,
#'     dimnames = list(statesNames, statesNames)))
#' impliedTimescales(mc)
#'
#' @exportMethod impliedTimescales
setGeneric("impliedTimescales", function(object) standardGeneric("impliedTimescales"))

#' @rdname impliedTimescales
setMethod("impliedTimescales", "markovchain", function(object) {
  moduli <- .nonTrivialEigenvalueModuli(object, "impliedTimescales")
  if (length(moduli) == 0L) {
    return(stats::setNames(numeric(0), character(0)))
  }

  tol <- sqrt(.Machine$double.eps)
  tau <- ifelse(
    moduli >= 1 - tol, Inf,
    ifelse(moduli <= tol, 0, -1 / log(moduli))
  )
  names(tau) <- paste0("tau", seq(2L, length(tau) + 1L))
  tau
})
