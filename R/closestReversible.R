#' Closest reversible approximation of a Markov chain
#'
#' Finds the Markov chain closest to a given one among those that are
#' reversible with respect to a fixed stationary distribution.
#'
#' For a row-stochastic transition matrix \eqn{P} with stationary
#' distribution \eqn{\pi}, define the \emph{time reversal} \eqn{P^*} by
#' \deqn{P^*_{ij} = \frac{\pi_j P_{ji}}{\pi_i}.}
#' \eqn{P^*} is the transition matrix of the same chain run backwards in
#' time, and \eqn{P} is reversible exactly when \eqn{P = P^*}. The
#' approximation returned is the \emph{additive reversibilization}
#' \deqn{R = \tfrac{1}{2}\left(P + P^*\right),}
#' which is stochastic, non-negative, and reversible with respect to the
#' same \eqn{\pi} (see Details).
#'
#' @param object A \code{markovchain} object representing a finite,
#'   irreducible discrete-time Markov chain.
#' @param stationaryDistribution Optional numeric vector giving the
#'   stationary distribution \eqn{\pi} to make the result reversible with
#'   respect to. It must have one strictly positive entry per state, sum to
#'   one (it is normalized if it does not), and be stationary for
#'   \code{object}, i.e. satisfy \eqn{\pi P = \pi}; otherwise an error is
#'   raised. The default, \code{NULL}, uses the chain's own unique
#'   stationary distribution from \code{\link{steadyStates}}.
#' @param tolerance A single finite non-negative number, used when checking
#'   that a supplied \code{stationaryDistribution} really is stationary. It
#'   is a numerical-noise tolerance, not a modelling one.
#'
#' @return A named list with four elements:
#'   \describe{
#'     \item{\code{chain}}{The approximating \code{markovchain} object
#'       \eqn{R}, with the same states and the same row/column-stochastic
#'       storage convention as \code{object}.}
#'     \item{\code{stationaryDistribution}}{The \eqn{\pi} used, as a named
#'       numeric vector.}
#'     \item{\code{distance}}{The distance \eqn{\|P-R\|_\pi} actually
#'       minimized (see Details).}
#'     \item{\code{frobeniusDistance}}{The plain Frobenius distance
#'       \eqn{\|P-R\|_F}, reported for convenience. It is \emph{not} the
#'       quantity being minimized.}
#'   }
#'
#' @details
#' \strong{In what sense is this the closest chain?} Work in the space
#' \eqn{\ell^2(\pi)} of functions on the states with inner product
#' \eqn{\langle f,g\rangle_\pi=\sum_i \pi_i f_i g_i}. Reversible chains are
#' exactly the self-adjoint operators on that space, and they form a linear
#' subspace. The matching Hilbert-Schmidt inner product on operators is
#' \deqn{\langle A,B\rangle_\pi = \sum_{i,j} \frac{\pi_i}{\pi_j} A_{ij}B_{ij},
#'   \qquad \|A\|_\pi^2 = \sum_{i,j} \frac{\pi_i}{\pi_j} A_{ij}^2,}
#' and \eqn{A \mapsto A^*} is an isometric involution for it. The orthogonal
#' projection onto the fixed points of such an involution is the average of
#' a point and its image, so \eqn{R=(P+P^*)/2} is the closest
#' \eqn{\pi}-reversible matrix to \eqn{P} in \eqn{\|\cdot\|_\pi}. The
#' stochastic and non-negativity constraints come for free: \eqn{P^*} has
#' row sums \eqn{\sum_j \pi_j P_{ji}/\pi_i = (\pi P)_i/\pi_i = 1} because
#' \eqn{\pi} is stationary, and both \eqn{P} and \eqn{P^*} are
#' non-negative, so the minimizer over the subspace already lies in the set
#' of transition matrices and no constrained optimization is needed.
#'
#' \strong{What this function does not do.} The minimization is over
#' reversible chains \emph{with \eqn{\pi} held fixed}, in the
#' \eqn{\pi}-weighted norm above. Two related problems are different and are
#' not solved here:
#' \itemize{
#'   \item Minimizing the \emph{plain} Frobenius distance
#'     \eqn{\|P-R\|_F} with \eqn{\pi} fixed. The involution
#'     \eqn{A\mapsto A^*} is not an isometry for that norm, so \eqn{R} is
#'     generally not its minimizer; \code{frobeniusDistance} is reported
#'     only as a descriptive figure.
#'   \item Letting the stationary distribution \emph{vary}, i.e. finding
#'     the reversible chain nearest to \eqn{P} over all choices of
#'     \eqn{\pi}. That is a genuinely harder constrained optimization
#'     problem, studied by Nielsen and Weber (2015), and it needs a
#'     numerical optimizer rather than a closed form. If you need it,
#'     supply candidate distributions through
#'     \code{stationaryDistribution} and compare \code{distance} values, or
#'     use a dedicated implementation.
#' }
#'
#' \strong{Properties worth knowing.} \eqn{R} has the same stationary
#' distribution \eqn{\pi} as \eqn{P}, and it preserves the support pattern
#' in the symmetrized sense: \eqn{R_{ij}>0} whenever \eqn{P_{ij}>0} or
#' \eqn{P_{ji}>0}. It may therefore allow transitions the original chain
#' forbids, which is inherent to making a chain reversible rather than a
#' defect of this construction. If \code{object} is already reversible,
#' \eqn{R=P} and \code{distance} is zero (up to rounding). Fill (1991)
#' introduces this construction, alongside the \emph{multiplicative}
#' reversibilization \eqn{PP^*}, which is a different object and is not
#' computed here.
#'
#' Only irreducibility is required, not aperiodicity. Irreducibility
#' guarantees both a unique \eqn{\pi} and \eqn{\pi_i>0} for every state, so
#' the division defining \eqn{P^*} is always safe.
#'
#' The implementation calls \code{\link{steadyStates}} at most once and is
#' then \eqn{O(n^2)} in time and memory for a dense \eqn{n}-state
#' transition matrix; no eigendecomposition or optimization is involved.
#'
#' @references
#' Fill, J. A. (1991). Eigenvalue bounds on convergence to stationarity for
#' nonreversible Markov chains, with an application to the exclusion
#' process. \emph{The Annals of Applied Probability}, 1(1).
#' \doi{10.1214/aoap/1177005981}
#'
#' Nielsen, A. and Weber, M. (2015). Computing the nearest reversible
#' Markov chain. \emph{Numerical Linear Algebra with Applications}, 22.
#' \doi{10.1002/nla.1967}
#'
#' @seealso \code{\link{is.reversible}}, \code{\link{steadyStates}},
#'   \code{\link{is.irreducible}}
#'
#' @examples
#' # A directed 3-cycle is as far from reversible as a chain gets: it only
#' # ever moves one way round. Its closest reversible approximation is the
#' # undirected random walk on the same triangle.
#' statesNames <- c("a", "b", "c")
#' cycle3 <- new("markovchain", states = statesNames,
#'   transitionMatrix = matrix(c(0, 1, 0,
#'                               0, 0, 1,
#'                               1, 0, 0), byrow = TRUE, nrow = 3,
#'                             dimnames = list(statesNames, statesNames)))
#' is.reversible(cycle3)
#'
#' approximation <- closestReversible(cycle3)
#' approximation$chain
#' is.reversible(approximation$chain)
#' approximation$distance
#'
#' @exportMethod closestReversible
setGeneric("closestReversible", function(object, stationaryDistribution = NULL,
                                         tolerance = sqrt(.Machine$double.eps)) {
  standardGeneric("closestReversible")
})

#' @rdname closestReversible
setMethod("closestReversible", "markovchain",
          function(object, stationaryDistribution = NULL,
                   tolerance = sqrt(.Machine$double.eps)) {
  if (length(tolerance) != 1L || !is.numeric(tolerance) ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be a single finite non-negative number.")
  }
  if (!is.irreducible(object)) {
    stop("closestReversible is defined here only for irreducible Markov chains.")
  }

  # Work row-stochastic internally; the result is put back into the input's
  # own orientation at the end.
  P <- as.matrix(object@transitionMatrix)
  if (!object@byrow) {
    P <- t(P)
  }
  n <- nrow(P)
  if (n != ncol(P) || any(!is.finite(P))) {
    stop("The transition matrix must be square and finite.")
  }

  stateNames <- states(object)

  if (is.null(stationaryDistribution)) {
    pi <- as.numeric(steadyStates(object))
    if (length(pi) != n || any(!is.finite(pi)) || sum(pi) <= 0) {
      stop("Unable to obtain a valid stationary distribution.")
    }
    pi[pi < 0] <- 0
    pi <- pi / sum(pi)
  } else {
    pi <- as.numeric(stationaryDistribution)
    if (length(pi) != n || any(!is.finite(pi))) {
      stop("stationaryDistribution must be a finite numeric vector with one entry per state.")
    }
    if (any(pi <= 0)) {
      stop(paste0(
        "stationaryDistribution must be strictly positive: the time reversal ",
        "divides by pi_i, which is undefined for a state of probability zero."
      ))
    }
    pi <- pi / sum(pi)
    # pi must be stationary for P, otherwise the time reversal below is not
    # even a transition matrix (its rows would not sum to one).
    drift <- max(abs(as.numeric(pi %*% P) - pi))
    if (drift > tolerance) {
      stop(paste0(
        "stationaryDistribution is not stationary for this chain: ",
        "max|pi P - pi| = ", format(drift, digits = 3),
        " exceeds tolerance. The closest reversible chain is only defined ",
        "with respect to a distribution the chain actually preserves."
      ))
    }
  }

  # Time reversal: Pstar[i, j] = pi[j] * P[j, i] / pi[i]. The sweep scales
  # column j of t(P) by pi[j]; dividing by pi then scales row i by 1/pi[i].
  Pstar <- sweep(t(P), 2L, pi, FUN = "*") / pi

  R <- (P + Pstar) / 2
  # Rows sum to one already; renormalize only to absorb rounding.
  R <- R / rowSums(R)
  dimnames(R) <- list(stateNames, stateNames)

  difference <- P - R
  # The pi-weighted Hilbert-Schmidt norm this construction minimizes:
  # ||A||^2 = sum_ij (pi_i / pi_j) A_ij^2. Scaling rows by pi_i and columns
  # by 1/pi_j gives that weight without forming an n-by-n weight matrix.
  weighted <- sweep(pi * difference^2, 2L, pi, FUN = "/")
  distance <- sqrt(sum(weighted))
  frobeniusDistance <- sqrt(sum(difference^2))

  if (!object@byrow) {
    R <- t(R)
  }

  list(
    chain = new("markovchain",
                states = stateNames,
                byrow = object@byrow,
                transitionMatrix = R,
                name = paste0(object@name, " (closest reversible)")),
    stationaryDistribution = stats::setNames(pi, stateNames),
    distance = distance,
    frobeniusDistance = frobeniusDistance
  )
})
