#' Closest reversible approximation of a Markov chain
#'
#' Finds the Markov chain closest to a given one among the reversible
#' chains, either with respect to a fixed stationary distribution or, with
#' \code{stationaryDistribution = "free"}, over all stationary
#' distributions.
#'
#' For a row-stochastic transition matrix \eqn{P} with stationary
#' distribution \eqn{\pi}, define the \emph{time reversal} \eqn{P^*} by
#' \deqn{P^*_{ij} = \frac{\pi_j P_{ji}}{\pi_i}.}
#' \eqn{P^*} is the transition matrix of the same chain run backwards in
#' time, and \eqn{P} is reversible exactly when \eqn{P = P^*}. With the
#' default \code{norm = "pi"} the approximation returned is the
#' \emph{additive reversibilization}
#' \deqn{R = \tfrac{1}{2}\left(P + P^*\right),}
#' which is stochastic, non-negative, and reversible with respect to the
#' same \eqn{\pi} (see Details). With \code{norm = "frobenius"} it is the
#' reversible chain nearest to \eqn{P} in the plain Frobenius norm, the
#' problem solved by Nielsen and Weber (2015).
#'
#' @param object A \code{markovchain} object representing a finite,
#'   irreducible discrete-time Markov chain.
#' @param stationaryDistribution Optional numeric vector giving the
#'   stationary distribution \eqn{\pi} to make the result reversible with
#'   respect to, with one strictly positive entry per state (it is
#'   normalized if it does not sum to one). With \code{norm = "pi"} it must
#'   be stationary for \code{object}, i.e. satisfy \eqn{\pi P = \pi},
#'   otherwise an error is raised; with \code{norm = "frobenius"} any
#'   positive distribution is allowed, and the result is the chain nearest to
#'   \eqn{P} among those reversible with respect to it. The default,
#'   \code{NULL}, uses the chain's own unique stationary distribution from
#'   \code{\link{steadyStates}}. With \code{norm = "frobenius"} the string
#'   \code{"free"} lets the stationary distribution vary as well (see
#'   Details).
#' @param tolerance A single finite non-negative number, used when checking
#'   that a supplied \code{stationaryDistribution} really is stationary
#'   (\code{norm = "pi"} only). It is a numerical-noise tolerance, not a
#'   modelling one.
#' @param norm The distance to minimize: \code{"pi"} (the default), the
#'   \eqn{\pi}-weighted Hilbert-Schmidt norm, whose minimizer is the
#'   additive reversibilization above; or \code{"frobenius"}, the plain
#'   Frobenius norm \eqn{\|P - R\|_F}, which needs the \pkg{quadprog}
#'   package.
#' @param nstart With \code{stationaryDistribution = "free"}, the number of
#'   starting points of the optimization (at least one); ignored otherwise.
#'   The first starts from the stationary distribution of \code{object} (or
#'   the uniform distribution if it is not unique), the second from the
#'   uniform distribution, the others from random distributions drawn with
#'   the session's random number generator.
#'
#' @return A named list with five elements:
#'   \describe{
#'     \item{\code{chain}}{The approximating \code{markovchain} object
#'       \eqn{R}, with the same states and the same row/column-stochastic
#'       storage convention as \code{object}.}
#'     \item{\code{stationaryDistribution}}{The \eqn{\pi} with respect to
#'       which \eqn{R} is reversible, as a named numeric vector: the one used
#'       or, with \code{"free"}, the one found.}
#'     \item{\code{distance}}{The distance actually minimized:
#'       \eqn{\|P-R\|_\pi} for \code{norm = "pi"}, \eqn{\|P-R\|_F} for
#'       \code{norm = "frobenius"} (see Details).}
#'     \item{\code{frobeniusDistance}}{The plain Frobenius distance
#'       \eqn{\|P-R\|_F}, reported for convenience; with
#'       \code{norm = "pi"} it is \emph{not} the quantity being minimized.}
#'     \item{\code{norm}}{The norm used.}
#'   }
#'   With \code{stationaryDistribution = "free"} a sixth element,
#'   \code{starts}, gives for each starting point the Frobenius distance
#'   reached.
#'
#' @details
#' \strong{In what sense is this the closest chain?} With the default
#' \code{norm = "pi"}, work in the space
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
#' \strong{The Frobenius norm with a fixed \eqn{\pi}.} With
#' \code{norm = "frobenius"} the result minimizes
#' \eqn{\|P - R\|_F^2 = \sum_{i,j}(P_{ij} - R_{ij})^2} over the transition
#' matrices that satisfy detailed balance
#' \eqn{\pi_i R_{ij} = \pi_j R_{ji}}. This is the problem of Nielsen and
#' Weber (2015). Since \eqn{A\mapsto A^*} is not an isometry of this norm,
#' the additive reversibilization is generally not its solution, and the
#' non-negativity constraints can be active. The problem is a convex
#' quadratic program: with the \eqn{n(n-1)/2} entries \eqn{R_{ij}}, \eqn{i<j},
#' as unknowns, \eqn{R_{ji} = \pi_i R_{ij}/\pi_j} and
#' \eqn{R_{ii} = 1 - \sum_{j\ne i}R_{ij}}, the objective is strictly convex
#' and the constraints are \eqn{R_{ij} \ge 0} and \eqn{R_{ii} \ge 0}. It
#' is solved exactly by the dual active-set method of
#' \code{quadprog::solve.QP}, in \eqn{O(n^6)} time, which is practical up to
#' about 50 states. \eqn{\pi} need not be stationary for \eqn{P}; the
#' result is reversible, hence stationary, with respect to it.
#'
#' \strong{Letting \eqn{\pi} vary.} With
#' \code{stationaryDistribution = "free"} the minimization is over all
#' reversible chains, whatever their stationary distribution. Every such
#' chain can be written as \eqn{R_{ij} = S_{ij}/\sum_k S_{ik}} with
#' \eqn{S} symmetric and non-negative (\eqn{S_{ij} = \pi_i R_{ij}}, the
#' probability flows), so the problem becomes a smooth minimization over
#' \eqn{S \ge 0}, solved with \code{optim(method = "L-BFGS-B")} and the
#' analytic gradient. The set of reversible chains is not convex, so the
#' result is a local minimum: \code{nstart} starting points are tried, and
#' each local solution is polished by solving the convex problem above at its
#' own \eqn{\pi}, which can only decrease the distance. Nielsen and Weber
#' (2015) keep \eqn{\pi} fixed; this extension is specific to the package.
#' The distance found is never larger than the one obtained with the
#' stationary distribution of \code{object}.
#'
#' \strong{What this function does not do.} With \code{norm = "pi"} the
#' minimization holds \eqn{\pi} fixed: the \eqn{\pi}-weighted norm itself
#' depends on \eqn{\pi}, so letting it vary is only offered for the
#' Frobenius norm. Also, the weighted Frobenius norm of PyDTMC's
#' \code{closest_reversible(weighted = TRUE)} is not offered. PyDTMC solves
#' the unweighted problem with an iterative solver, whose solutions agree with
#' the ones here to about \eqn{10^{-5}} and have a slightly larger distance.
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
#' # The nearest reversible chain in the Frobenius norm (needs quadprog),
#' # first with the stationary distribution held fixed, then free.
#' if (requireNamespace("quadprog", quietly = TRUE)) {
#'   P <- new("markovchain", states = letters[1:4], transitionMatrix = matrix(
#'     c(0.1, 0.6, 0.2, 0.1, 0.05, 0.15, 0.7, 0.1,
#'       0.3, 0.05, 0.25, 0.4, 0.5, 0.2, 0.1, 0.2), 4, byrow = TRUE,
#'     dimnames = list(letters[1:4], letters[1:4])))
#'   closestReversible(P)$frobeniusDistance
#'   closestReversible(P, norm = "frobenius")$distance
#'   set.seed(1)
#'   free <- closestReversible(P, "free", norm = "frobenius")
#'   free$distance
#'   free$stationaryDistribution
#' }
#'
#' @exportMethod closestReversible
# Closest reversible chain to P: the one that satisfies detailed balance pi_i r_ij = pi_j r_ji.
# norm = "pi": the average of P and its time reversal (closed form).
# norm = "frobenius": the quadratic program of Nielsen and Weber (2015), with pi fixed or, with
# stationaryDistribution = "free", optimized as well.
setGeneric("closestReversible", function(object, stationaryDistribution = NULL,
                                         tolerance = sqrt(.Machine$double.eps),
                                         norm = c("pi", "frobenius"), nstart = 5L) {
  standardGeneric("closestReversible")
})

# Internal helper: nearest matrix to P (row-stochastic, n x n) in the Frobenius
# norm among the transition matrices reversible with respect to pi (> 0).
# Unknowns: y_k = R[i, j] for the pairs i < j. Detailed balance fixes
# R[j, i] = rho_k y_k with rho_k = pi_i / pi_j, and R[i, i] = 1 - (row sum of
# the off-diagonal entries) = 1 - C[i, ] %*% y, where C[i, k] is 1 if i is the
# first state of pair k and rho_k if it is the second. The objective
#   sum_k (y_k - P_ij)^2 + (rho_k y_k - P_ji)^2 + sum_i (1 - P_ii - C[i, ] y)^2
# is y' A y - 2 b' y + constant, with A = diag(1 + rho^2) + C'C positive
# definite, and the constraints are y >= 0 and C y <= 1.
.nearestReversibleFixedPi <- function(P, pi) {
  n <- nrow(P)
  if (n < 2L) {
    return(P) # a single state is always reversible
  }
  pairs <- which(upper.tri(P), arr.ind = TRUE)
  i <- pairs[, 1L]
  j <- pairs[, 2L]
  m <- length(i)
  rho <- pi[i] / pi[j]
  C <- matrix(0, n, m)
  C[cbind(i, seq_len(m))] <- 1
  C[cbind(j, seq_len(m))] <- rho
  A <- diag(1 + rho^2, m) + crossprod(C)
  b <- P[cbind(i, j)] + rho * P[cbind(j, i)] + as.numeric(crossprod(C, 1 - diag(P)))
  # solve.QP minimizes -d'y + y'Dy/2 subject to t(Amat) y >= bvec
  sol <- quadprog::solve.QP(Dmat = 2 * A, dvec = 2 * b,
                            Amat = cbind(diag(m), -t(C)),
                            bvec = c(rep(0, m), rep(-1, n)))
  y <- pmax(sol$solution, 0) # remove rounding below zero
  R <- matrix(0, n, n)
  R[cbind(i, j)] <- y
  R[cbind(j, i)] <- rho * y
  diag(R) <- pmax(1 - rowSums(R), 0)
  R / rowSums(R)
}

# Internal helper: nearest reversible matrix to P in the Frobenius norm with
# the stationary distribution free. A reversible chain is R = S / rowSums(S)
# with S symmetric and non-negative (S = diag(pi) R, the flows). The unknowns
# are the upper triangle of S, diagonal included; the diagonal has a tiny
# positive lower bound so that no row of S can vanish. F(S) = ||S / r - P||^2
# is homogeneous of degree zero in S, and its gradient, with E = R - P and
# r = rowSums(S), is dF / dS[i, j] = 2 (E[i, j] - sum_k E[i, k] R[i, k]) / r_i
# for an unconstrained entry, summed over (i, j) and (j, i) for the symmetric
# unknowns.
.nearestReversibleFreePi <- function(P, startPis) {
  n <- nrow(P)
  upper <- upper.tri(P, diag = TRUE)               # positions of the free unknowns u
  onDiagonal <- (row(P) == col(P))[upper]
  # build(u): the symmetric matrix S with upper triangle u, and the reversible chain
  # R = S / rowSums(S) it defines (r = rowSums(S) is proportional to its stationary law).
  build <- function(u) {
    S <- matrix(0, n, n)
    S[upper] <- u
    S <- S + t(S) - diag(diag(S), n)               # mirror below the diagonal
    r <- rowSums(S)
    list(R = S / r, r = r)
  }
  objective <- function(u) sum((build(u)$R - P)^2)   # squared Frobenius distance to P
  # analytic gradient of the objective with respect to u (formula in the header above)
  gradient <- function(u) {
    b <- build(u)
    E <- b$R - P
    H <- 2 * (E - rowSums(E * b$R)) / b$r
    G <- H + t(H)                                  # an off-diagonal u enters S twice
    diag(G) <- diag(H)
    G[upper]
  }
  lower <- ifelse(onDiagonal, 1e-12, 0)            # bounds: S >= 0, diagonal slightly positive
  # Local optimization from each starting stationary distribution; the best result is kept.
  results <- lapply(startPis, function(pi0) {
    R0 <- .nearestReversibleFixedPi(P, pi0)        # exact solution for the starting pi
    S0 <- pi0 * R0                                 # its flow matrix, the starting point
    opt <- stats::optim(S0[upper], objective, gradient, method = "L-BFGS-B",
                        lower = lower, control = list(maxit = 10000L, factr = 10))
    found <- build(opt$par)
    pi1 <- found$r / sum(found$r)
    # polish: the exact solution of the convex problem at the pi just found
    # is never worse than the local solution, which is feasible for it
    R1 <- .nearestReversibleFixedPi(P, pi1)
    if (sum((R1 - P)^2) > sum((found$R - P)^2)) {
      R1 <- found$R
    }
    list(R = R1, pi = pi1, distance = sqrt(sum((R1 - P)^2)), convergence = opt$convergence)
  })
  best <- which.min(vapply(results, `[[`, numeric(1), "distance"))
  out <- results[[best]]
  out$starts <- data.frame(start = seq_along(results),
                           distance = vapply(results, `[[`, numeric(1), "distance"),
                           convergence = vapply(results, `[[`, integer(1), "convergence"))
  out
}

#' @rdname closestReversible
setMethod("closestReversible", "markovchain",
          function(object, stationaryDistribution = NULL,
                   tolerance = sqrt(.Machine$double.eps),
                   norm = c("pi", "frobenius"), nstart = 5L) {
  # Step 1: validate the arguments and the chain.
  norm <- match.arg(norm)
  if (length(tolerance) != 1L || !is.numeric(tolerance) ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be a single finite non-negative number.")
  }
  free <- is.character(stationaryDistribution)
  if (free) {
    if (!identical(stationaryDistribution, "free")) {
      stop("stationaryDistribution must be NULL, a numeric vector or \"free\".")
    }
    if (norm != "frobenius") {
      stop("stationaryDistribution = \"free\" requires norm = \"frobenius\": the ",
           "pi-weighted norm itself depends on the stationary distribution.")
    }
    if (length(nstart) != 1L || !is.numeric(nstart) || is.na(nstart) ||
        nstart != as.integer(nstart) || nstart < 1L) {
      stop("nstart must be a single integer of at least 1.")
    }
  }
  if (norm == "frobenius" && !requireNamespace("quadprog", quietly = TRUE)) {
    stop("norm = \"frobenius\" requires the 'quadprog' package. ",
         "Install it with install.packages(\"quadprog\").")
  }
  # The chain's own stationary distribution is needed for norm = "pi" and for
  # the default of norm = "frobenius"; only then must it be unique.
  if ((norm == "pi" || is.null(stationaryDistribution)) && !is.irreducible(object)) {
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
  starts <- NULL

  # Step 2: obtain pi, computed from the chain, supplied by the user or (free) optimized.
  if (free) {
    pi0 <- if (is.irreducible(object)) {
      as.numeric(steadyStates(object))
    } else {
      rep(1 / n, n)
    }
    pi0 <- pmax(pi0, 0)
    pi0 <- pi0 / sum(pi0)
    # The problem is not convex in pi, so several starting points are used: the chain's own
    # stationary law, the uniform law, then random laws (uniform on the simplex) up to nstart.
    startPis <- c(list(pi0, rep(1 / n, n)),
                  lapply(seq_len(max(0L, as.integer(nstart) - 2L)), function(k) {
                    g <- stats::rexp(n)
                    g / sum(g)
                  }))[seq_len(nstart)]
    found <- .nearestReversibleFreePi(P, startPis)
    R <- found$R
    pi <- found$pi
    starts <- found$starts
  } else {
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
      # With norm = "pi", pi must be stationary for P, otherwise the time
      # reversal below is not even a transition matrix (its rows would not
      # sum to one). The Frobenius problem accepts any positive pi.
      drift <- max(abs(as.numeric(pi %*% P) - pi))
      if (norm == "pi" && drift > tolerance) {
        stop(paste0(
          "stationaryDistribution is not stationary for this chain: ",
          "max|pi P - pi| = ", format(drift, digits = 3),
          " exceeds tolerance. The closest reversible chain is only defined ",
          "with respect to a distribution the chain actually preserves."
        ))
      }
    }

    if (norm == "pi") {
      # Step 3: time reversal P*[i,j] = pi_j P[j,i] / pi_i.
      # Time reversal: Pstar[i, j] = pi[j] * P[j, i] / pi[i]. The sweep scales
      # column j of t(P) by pi[j]; dividing by pi then scales row i by 1/pi[i].
      Pstar <- sweep(t(P), 2L, pi, FUN = "*") / pi

      # Step 4: the average of P and P* is reversible and is the closest in the pi-weighted norm.
      R <- (P + Pstar) / 2
      # Rows sum to one already; renormalize only to absorb rounding.
      R <- R / rowSums(R)
    } else {
      # Steps 3-4 (Frobenius): the convex quadratic program of Nielsen and Weber.
      R <- .nearestReversibleFixedPi(P, pi)
    }
  }
  dimnames(R) <- list(stateNames, stateNames)

  # Step 5: report how far the result is from the original chain (two different norms).
  difference <- P - R
  frobeniusDistance <- sqrt(sum(difference^2))
  if (norm == "pi") {
    # The pi-weighted Hilbert-Schmidt norm this construction minimizes:
    # ||A||^2 = sum_ij (pi_i / pi_j) A_ij^2. Scaling rows by pi_i and columns
    # by 1/pi_j gives that weight without forming an n-by-n weight matrix.
    weighted <- sweep(pi * difference^2, 2L, pi, FUN = "/")
    distance <- sqrt(sum(weighted))
  } else {
    distance <- frobeniusDistance
  }

  # Return in the same orientation (rows or columns) as the input object.
  if (!object@byrow) {
    R <- t(R)
  }

  out <- list(
    chain = new("markovchain",
                states = stateNames,
                byrow = object@byrow,
                transitionMatrix = R,
                name = paste0(object@name, " (closest reversible)")),
    stationaryDistribution = stats::setNames(pi, stateNames),
    distance = distance,
    frobeniusDistance = frobeniusDistance,
    norm = norm
  )
  if (!is.null(starts)) {
    out$starts <- starts
  }
  out
})
