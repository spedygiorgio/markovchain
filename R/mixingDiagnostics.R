# These two functions both ask "how does this chain approach its
# stationary distribution?", but they answer it in different ways and
# therefore require different assumptions:
#
#   * is.reversible() only checks an algebraic property (detailed balance)
#     that makes sense for ANY irreducible chain, periodic or not.
#   * mixingTime() asks how many steps are needed before the chain's
#     distribution is numerically close to stationary. For a periodic
#     chain this question has no good answer -- the distribution never
#     settles down, it keeps oscillating -- so mixingTime() additionally
#     requires aperiodicity (see its own documentation for why).

#' Check whether a Markov chain is reversible
#'
#' Checks whether a finite, irreducible discrete-time Markov chain is
#' reversible with respect to its (unique) stationary distribution, i.e.
#' whether it satisfies the detailed balance equations.
#'
#' A chain with transition matrix \eqn{P} and stationary distribution
#' \eqn{\pi} is reversible if
#' \deqn{\pi_i P_{ij} = \pi_j P_{ji} \quad \text{for every } i,j.}
#' Intuitively, if you started the chain from \eqn{\pi} and watched a long
#' run of it, running the recorded sequence of states backwards would look
#' statistically identical to running it forwards: at stationarity, the
#' "flow" of probability from \eqn{i} to \eqn{j} exactly balances the flow
#' from \eqn{j} to \eqn{i}.
#'
#' @param object A \code{markovchain} object representing a finite,
#'   irreducible discrete-time Markov chain.
#' @param tolerance A single finite non-negative number. Detailed balance is
#'   accepted as holding when every pair \eqn{(i,j)} satisfies
#'   \eqn{|\pi_i P_{ij} - \pi_j P_{ji}| \le \code{tolerance}}. The default is
#'   a small numerical-noise tolerance, not a modelling tolerance: it exists
#'   to absorb floating-point rounding in the eigendecomposition-based
#'   \code{\link{steadyStates}} computation, not to declare "almost
#'   reversible" chains reversible.
#'
#' @return A single logical value, \code{TRUE} or \code{FALSE}.
#'
#' @details
#' Only irreducibility is required, not aperiodicity: detailed balance is a
#' purely algebraic condition on \eqn{P} and \eqn{\pi} and is perfectly well
#' defined for periodic chains too. For example, a simple random walk on
#' any undirected graph (moving to a uniformly random neighbour) is always
#' reversible, whether or not it happens to be periodic.
#'
#' A 2-state irreducible chain is always reversible: with only two states,
#' the single detailed balance equation \eqn{\pi_1 P_{12} = \pi_2 P_{21}} is
#' just a restatement of the stationarity equation \eqn{\pi P = \pi}, so it
#' holds automatically.
#'
#' Every reversible chain has a real spectrum (all eigenvalues of \eqn{P}
#' are real), which is why \code{\link{slem}} and \code{\link{spectralGap}}
#' are especially easy to interpret for reversible chains: there are no
#' complex-conjugate eigenvalue pairs to reason about.
#'
#' The implementation calls \code{\link{steadyStates}} once, then compares
#' the two triangles of the flow matrix \eqn{\pi_i P_{ij}}. Its time
#' complexity is dominated by \code{steadyStates()}, plus an additional
#' \eqn{O(n^2)} comparison for a dense \eqn{n}-state transition matrix. It
#' supports both row- and column-stochastic storage.
#'
#' @references
#' Norris, J. R. (1998). \emph{Markov Chains}. Cambridge University Press.
#'
#' Levin, D. A. and Peres, Y. (2017). \emph{Markov Chains and Mixing Times},
#' 2nd edition. American Mathematical Society.
#'
#' @seealso \code{\link{steadyStates}}, \code{\link{is.irreducible}},
#'   \code{\link{slem}}, \code{\link{mixingTime}}
#'
#' @examples
#' # A random walk on a triangle is reversible.
#' statesNames <- c("a", "b", "c")
#' triangle <- new("markovchain", states = statesNames,
#'   transitionMatrix = matrix(c(0, 0.5, 0.5,
#'                               0.5, 0, 0.5,
#'                               0.5, 0.5, 0), byrow = TRUE, nrow = 3,
#'                             dimnames = list(statesNames, statesNames)))
#' is.reversible(triangle)
#'
#' # A directed cycle (states only move "forward") is not reversible: there
#' # is a net clockwise flow of probability at stationarity.
#' cycle3 <- new("markovchain", states = statesNames,
#'   transitionMatrix = matrix(c(0, 1, 0,
#'                               0, 0, 1,
#'                               1, 0, 0), byrow = TRUE, nrow = 3,
#'                             dimnames = list(statesNames, statesNames)))
#' is.reversible(cycle3)
#'
#' @exportMethod is.reversible
setGeneric("is.reversible", function(object, tolerance = sqrt(.Machine$double.eps)) {
  standardGeneric("is.reversible")
})

#' @rdname is.reversible
setMethod("is.reversible", "markovchain", function(object, tolerance = sqrt(.Machine$double.eps)) {
  if (length(tolerance) != 1L || !is.numeric(tolerance) ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be a single finite non-negative number.")
  }
  if (!is.irreducible(object)) {
    stop("is.reversible is defined here only for irreducible Markov chains.")
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
    return(TRUE) # a single-state chain trivially satisfies detailed balance
  }

  pi <- as.numeric(steadyStates(object))
  if (length(pi) != n || any(!is.finite(pi)) ||
      any(pi < -sqrt(.Machine$double.eps)) || sum(pi) <= 0) {
    stop("Unable to obtain a valid stationary distribution.")
  }
  pi[pi < 0] <- 0
  pi <- pi / sum(pi)

  # flow[i, j] = pi_i * P_ij. Multiplying the matrix by the vector `pi`
  # recycles pi down each column, i.e. row i is scaled by pi[i] -- exactly
  # the row-wise scaling detailed balance needs, with no explicit
  # n-by-n outer product ever constructed.
  flow <- pi * P

  max(abs(flow - t(flow))) <= tolerance
})

#' Mixing time of a Markov chain
#'
#' Estimates the total-variation mixing time of a finite, irreducible,
#' aperiodic (i.e. ergodic) discrete-time Markov chain: the smallest number
#' of steps after which the chain's distribution is within \code{epsilon} of
#' its stationary distribution, from every possible starting state.
#'
#' For a row-stochastic transition matrix \eqn{P} with stationary
#' distribution \eqn{\pi}, define the worst-case total variation distance
#' after \eqn{t} steps as
#' \deqn{d(t) = \max_i \tfrac{1}{2}\sum_j |P^t_{ij} - \pi_j|.}
#' The mixing time returned is
#' \deqn{t_{\mathrm{mix}}(\varepsilon) = \min\{t \ge 1 : d(t) \le \varepsilon\}.}
#'
#' @param object A \code{markovchain} object representing a finite,
#'   irreducible, aperiodic discrete-time Markov chain.
#' @param epsilon A single number strictly between \code{0} and \code{1}:
#'   the total-variation threshold that counts as "mixed". The classical
#'   default \code{0.25} follows Levin and Peres (2017); it is a
#'   conventional choice, not a universal constant.
#' @param maxIter A single positive integer: the largest \eqn{t} that will
#'   be tried before giving up. This is a safety limit, not a modelling
#'   parameter: it exists so that a chain which (numerically) mixes only
#'   extremely slowly reports a clear error instead of looping for an
#'   unbounded number of iterations.
#'
#' @return A single positive integer, the estimated mixing time
#'   \eqn{t_{\mathrm{mix}}(\varepsilon)}. For the trivial one-state chain,
#'   \code{0} is returned (it is its own stationary distribution).
#'
#' @details
#' Unlike \code{\link{slem}}, \code{\link{spectralGap}} and
#' \code{\link{impliedTimescales}}, \code{mixingTime()} requires
#' aperiodicity in addition to irreducibility. This is not an arbitrary
#' restriction carried over from another implementation: for a periodic
#' chain, \eqn{P^t(i, \cdot)} never converges to \eqn{\pi} at all (it keeps
#' cycling through a fixed set of distributions), so \eqn{d(t)} does not go
#' to zero and "the number of steps until \eqn{d(t) \le \varepsilon}" is
#' simply undefined for small enough \eqn{\varepsilon}. \code{slem()} and
#' \code{spectralGap()} remain meaningful for periodic chains because they
#' summarise the transition matrix's spectrum directly, without reference
#' to a limit that may not exist.
#'
#' The implementation repeatedly forms \eqn{P^{t+1} = P^t P} and checks
#' \eqn{d(t)} after each multiplication, starting from \eqn{t=1}, until the
#' threshold is met or \code{maxIter} is reached. Its time complexity is
#' \eqn{O(t_{\mathrm{mix}} \cdot n^3)} and its memory use is \eqn{O(n^2)}
#' for a dense \eqn{n}-state transition matrix: this is a direct,
#' easy-to-audit computation, not an asymptotically optimal one (a
#' repeated-squaring scheme would reach a single large power of \eqn{P}
#' faster, but would not let every intermediate \eqn{t} be checked against
#' \code{epsilon} along the way). It supports both row- and
#' column-stochastic storage.
#'
#' @references
#' Levin, D. A. and Peres, Y. (2017). \emph{Markov Chains and Mixing Times},
#' 2nd edition. American Mathematical Society.
#'
#' @seealso \code{\link{slem}}, \code{\link{spectralGap}},
#'   \code{\link{impliedTimescales}}, \code{\link{period}},
#'   \code{\link{is.irreducible}}
#'
#' @examples
#' statesNames <- c("a", "b")
#' mc <- new("markovchain",
#'   states = statesNames,
#'   transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
#'     byrow = TRUE, nrow = 2,
#'     dimnames = list(statesNames, statesNames)))
#' mixingTime(mc)
#' mixingTime(mc, epsilon = 0.01) # a tighter threshold needs more steps
#'
#' @exportMethod mixingTime
setGeneric("mixingTime", function(object, epsilon = 0.25, maxIter = 10000L) {
  standardGeneric("mixingTime")
})

#' @rdname mixingTime
setMethod("mixingTime", "markovchain", function(object, epsilon = 0.25, maxIter = 10000L) {
  if (length(epsilon) != 1L || !is.numeric(epsilon) || !is.finite(epsilon) ||
      epsilon <= 0 || epsilon >= 1) {
    stop("epsilon must be a single finite number strictly between 0 and 1.")
  }
  if (length(maxIter) != 1L || !is.numeric(maxIter) || !is.finite(maxIter) ||
      maxIter != as.integer(maxIter) || maxIter < 1L) {
    stop("maxIter must be a single positive integer.")
  }
  maxIter <- as.integer(maxIter)

  if (!is.irreducible(object)) {
    stop("mixingTime is defined here only for irreducible Markov chains.")
  }
  if (period(object) != 1L) {
    stop(paste0(
      "mixingTime is defined here only for aperiodic (ergodic) Markov ",
      "chains, because a periodic chain never converges to its ",
      "stationary distribution. See slem() or spectralGap() for a ",
      "spectral diagnostic that remains meaningful for periodic chains."
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
    return(0L)
  }

  pi <- as.numeric(steadyStates(object))
  if (length(pi) != n || any(!is.finite(pi)) ||
      any(pi < -sqrt(.Machine$double.eps)) || sum(pi) <= 0) {
    stop("Unable to obtain a valid stationary distribution.")
  }
  pi[pi < 0] <- 0
  pi <- pi / sum(pi)

  Pt <- P
  t <- 1L
  repeat {
    # Total variation distance from each row of P^t to pi, halved sum of
    # absolute differences; take the worst (largest) over starting states.
    distances <- 0.5 * rowSums(abs(sweep(Pt, 2L, pi, FUN = "-")))
    if (max(distances) <= epsilon) {
      return(t)
    }
    if (t >= maxIter) {
      stop(paste0(
        "mixingTime did not reach the epsilon threshold within maxIter = ",
        maxIter, " steps; the chain may mix very slowly, or maxIter may ",
        "need to be increased."
      ))
    }
    Pt <- Pt %*% P
    t <- t + 1L
  }
})
