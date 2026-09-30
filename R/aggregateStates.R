# Internal helpers for aggregateStates() ------------------------------------
#
# These implement the spectral-theoretic Kullback-Leibler (KL) state-space
# aggregation of Deng, Mehta and Meyn (2011): given a finite, irreducible,
# aperiodic Markov chain with transition matrix P and stationary distribution
# pi, and a 0/1 indicator matrix phi (n rows, one column per macro-state,
# each row summing to one) describing a partition of the n micro-states into
# macro-states, the induced macro-chain that best approximates P in the
# sense of minimizing the KL divergence rate between the two processes has
# transition matrix
#   Q = D^{-1} phi' Pi P phi,   D = diag(phi' Pi phi 1)
# where Pi = diag(pi). Two greedy heuristics build a good phi for a target
# number of macro-states k: "bottom-up" grows the partition one split at a
# time starting from a single macro-state (suited to a large reduction, i.e.
# k much smaller than n), and "top-down" repeatedly merges the pair of
# macro-states that increases the KL divergence least, starting from every
# micro-state on its own (suited to a small reduction, i.e. k close to n).
# Both use the same building block for choosing where to split/merge: the
# additive reversibilization (Fill, 1991; the same construction used by
# closestReversible()) of the relevant sub-chain, whose second eigenvector's
# sign pattern gives a near-optimal two-way cut (this is the standard
# spectral bisection heuristic, not a KL-optimal cut in itself; the KL
# divergence is only evaluated, and used to choose among candidates, after a
# cut is proposed).

# phi: n x z indicator matrix (z macro-states). Returns the z x z
# macro-transition matrix Q implied by phi, using the closed-form projection
# above. Shared by both the bottom-up and top-down heuristics below.
.aggregationQ <- function(P, pi, phi) {
  piDiag <- diag(pi, nrow = length(pi))
  z <- ncol(phi)
  qNum <- t(phi) %*% piDiag %*% P %*% phi
  qDen <- numeric(z)
  for (zi in seq_len(z)) {
    col <- phi[, zi]
    qDen[zi] <- as.numeric(t(col) %*% piDiag %*% col)
  }
  qNum / matrix(qDen, nrow = z, ncol = z)
}

# The information loss of representing (P, pi) by the aggregated (Q, phi),
# Deng et al.'s equation (8): the KL divergence rate between the original
# chain and the "lifted" macro-chain, minus the KL divergence rate already
# incurred just by aggregating the stationary distribution itself. Convention
# 0*log(0/y) = 0 (an all-zero row of P never contributes, including against
# a superQ/superPi of 0) is applied directly rather than via a NaN patch-up.
.klDivergenceAggregate <- function(P, pi, phi, Q) {
  superQ <- phi %*% Q %*% t(phi)
  ratioP <- ifelse(P == 0, 1, P / superQ)
  theta <- rowSums(P * log2(ratioP))
  tVal <- sum(pi * theta)

  superPi <- as.numeric(pi %*% phi %*% t(phi))
  psi <- log2(pi / superPi)
  uVal <- sum(pi * psi)

  tVal - uVal
}

# Additive reversibilization of M restricted to (and pi-weighted within) the
# states in idx, split by the sign of its second eigenvector. Returns a
# logical vector aligned with idx (TRUE/FALSE = the two sides of the cut).
# Falls back to a rank-based split (guaranteed non-empty on both sides) in
# the rare degenerate case where the eigenvector does not change sign.
.spectralSplit <- function(M, pi, idx) {
  Msub <- M[idx, idx, drop = FALSE]
  piSub <- diag(pi[idx], nrow = length(idx))
  AR <- 0.5 * (Msub + solve(piSub, t(Msub)) %*% piSub)
  eig <- eigen(AR)
  ord <- order(Mod(eig$values))
  secondIdx <- ord[length(ord) - 1L]
  evector <- Re(eig$vectors[, secondIdx])

  pos <- evector >= 0
  if (all(pos) || !any(pos)) {
    rankOrd <- order(evector)
    pos <- rep(FALSE, length(idx))
    pos[rankOrd[seq_len(length(idx) %/% 2L)]] <- TRUE
  }
  pos
}

# Bottom-up: split macro-state `index` of the current n x (ncol(phi))
# partition `phi` into two, returning the resulting n x (ncol(phi)+1)
# candidate partition, or NULL if that macro-state is a singleton (nothing
# to split).
.createBipartitionCandidate <- function(P, pi, phi, index) {
  v <- phi[, index]
  if (sum(v) <= 1) {
    return(NULL)
  }
  idx <- which(v > 0)
  pos <- .spectralSplit(P, pi, idx)

  v1 <- v; v1[idx] <- as.numeric(pos)
  v2 <- v; v2[idx] <- as.numeric(!pos)

  before <- if (index > 1L) phi[, seq_len(index - 1L), drop = FALSE] else NULL
  after <- if (index < ncol(phi)) phi[, seq.int(index + 1L, ncol(phi)), drop = FALSE] else NULL
  cbind(before, v1, v2, after)
}

# Bottom-up heuristic (Deng et al., Algorithm 1): start from the trivial
# one-macro-state partition and, at each step, split whichever macro-state
# yields the smallest resulting KL divergence, until k macro-states are
# reached.
.aggregateSpectralBottomUp <- function(P, pi, k) {
  n <- nrow(P)
  phi <- matrix(1, nrow = n, ncol = 1)
  Q <- matrix(1 / n, n, n)
  cur <- 1L

  while (cur < k) {
    bestPhi <- NULL
    bestQ <- NULL
    bestR <- Inf

    for (i in seq_len(ncol(phi))) {
      candidate <- .createBipartitionCandidate(P, pi, phi, i)
      if (is.null(candidate)) {
        next
      }
      candidateQ <- .aggregationQ(P, pi, candidate)
      candidateR <- .klDivergenceAggregate(P, pi, candidate, candidateQ)
      if (candidateR < bestR) {
        bestPhi <- candidate
        bestQ <- candidateQ
        bestR <- candidateR
      }
    }

    if (is.null(bestPhi)) {
      stop("Unable to split any further macro-state (all are singletons); try a smaller k or method = \"spectral-top-down\".")
    }
    phi <- bestPhi
    Q <- bestQ
    cur <- cur + 1L
  }

  list(Q = Q, phi = phi)
}

# The stationary distribution of a (possibly periodic) row-stochastic matrix
# Q, by simple Cesaro-averaged power iteration; used only to weight the
# spectral cuts proposed while growing the top-down candidate pairs below,
# never as the divergence's own pi (which is always the original chain's).
.calculateInvariant <- function(Q) {
  n <- nrow(Q)
  kappa <- rep(1 / n, n)
  theta <- as.numeric(kappa %*% Q)
  z <- 0L
  while (max(abs(kappa - theta)) > 1e-8 && z < 1000L) {
    kappa <- (kappa + theta) / 2
    theta <- as.numeric(kappa %*% Q)
    z <- z + 1L
  }
  theta
}

# Recursively spectral-bisect the current macro-state indices 1:ncol(Q)
# until every piece has at most two members, and return the size-2 pieces as
# candidate pairs to merge. A piece that splits down to a lone singleton
# (only possible when the piece being split has odd size) is simply not a
# usable merge candidate and is dropped -- there is nothing to pair it with
# yet at this recursion level, not a state lost from the chain: every state
# still belongs to some ancestor piece, and remains available the next time
# aggregateStates() rebuilds the candidate list from the (unchanged) states
# outside whichever pair actually gets merged.
.buildPairCandidates <- function(Q, pi) {
  n <- ncol(Q)
  worklist <- list(seq_len(n))
  pairs <- list()

  while (length(worklist) > 0L) {
    grp <- worklist[[1L]]
    worklist <- worklist[-1L]

    if (length(grp) <= 1L) {
      next
    }
    if (length(grp) == 2L) {
      pairs[[length(pairs) + 1L]] <- sort(grp)
      next
    }

    pos <- .spectralSplit(Q, pi, grp)
    worklist <- c(worklist, list(grp[pos]), list(grp[!pos]))
  }

  pairs
}

# Merge columns vi0 and vi1 (vi0 < vi1) of the current n x m assignment
# matrix `eta` into one, preserving every other column (and the relative
# order of the states within it) exactly. Returns the resulting n x (m-1)
# matrix together with the macro-chain it implies.
.calculateQTopDown <- function(P, pi, eta, vi0, vi1) {
  mergedCol <- pmax(eta[, vi0], eta[, vi1])
  keep <- setdiff(seq_len(ncol(eta)), c(vi0, vi1))
  before <- keep[keep < vi0]
  middle <- keep[keep > vi0 & keep < vi1]
  after <- keep[keep > vi1]

  phiNew <- cbind(eta[, before, drop = FALSE], mergedCol,
                  eta[, middle, drop = FALSE], eta[, after, drop = FALSE])

  list(Q = .aggregationQ(P, pi, phiNew), phi = phiNew)
}

# Top-down heuristic (Deng et al., Algorithm 2): start from every
# micro-state as its own macro-state and, at each step, merge whichever pair
# of current macro-states yields the smallest resulting KL divergence, until
# only k macro-states remain.
.aggregateSpectralTopDown <- function(P, pi, k) {
  n <- nrow(P)
  Q <- P
  eta <- diag(n)
  iterations <- n - k

  for (i in seq_len(iterations)) {
    qPi <- if (i == 1L) pi else .calculateInvariant(Q)
    pairs <- .buildPairCandidates(Q, qPi)

    divergences <- numeric(length(pairs))
    candidates <- vector("list", length(pairs))
    for (j in seq_along(pairs)) {
      res <- .calculateQTopDown(P, pi, eta, pairs[[j]][1L], pairs[[j]][2L])
      candidates[[j]] <- res
      divergences[j] <- .klDivergenceAggregate(P, pi, res$phi, res$Q)
    }

    best <- which.min(divergences)
    Q <- candidates[[best]]$Q
    eta <- candidates[[best]]$phi
  }

  list(Q = Q, phi = eta)
}

# Number of macro-states suggested by the eigengap heuristic: the number of
# eigenvalues (by modulus, including the trivial one) before the largest
# relative drop, restricted to the valid aggregateStates() range [2, n-1].
.chooseKByEigengap <- function(P) {
  n <- nrow(P)
  ev <- sort(Mod(eigen(P, only.values = TRUE)$values), decreasing = TRUE)
  gaps <- ev[-n] - ev[-1L]
  candidates <- seq(2L, n - 1L)
  candidates[which.max(gaps[candidates])]
}

#' Aggregate a Markov chain's state space by Kullback-Leibler minimization
#'
#' Reduces the state space of a finite, irreducible, aperiodic Markov chain
#' to \code{k} macro-states by the spectral-theoretic method of Deng, Mehta
#' and Meyn (2011): the macro-chain returned is the one whose "lifted"
#' behavior (each macro-state visit standing in for its micro-states,
#' weighted by their share of the stationary distribution) is closest, in
#' Kullback-Leibler divergence rate, to the original chain.
#'
#' @param object A \code{markovchain} object representing a finite,
#'   irreducible, aperiodic discrete-time Markov chain with at least 3
#'   states.
#' @param k The number of macro-states to reduce to, an integer between 2
#'   and the number of states minus 1. The default, \code{NULL}, chooses
#'   \code{k} automatically via the eigengap heuristic: the transition
#'   matrix's eigenvalues are sorted by modulus, and \code{k} is set to the
#'   number of eigenvalues before the largest relative drop -- a standard,
#'   parameter-free way to guess how many "slow", well-separated modes the
#'   chain has. This is a genuine automatic selection, unlike PyDTMC's
#'   \code{adaptive}, which only picks which of the two algorithms below to
#'   run for a \code{k} the caller must still supply.
#' @param method One of \code{"adaptive"} (the default), \code{"spectral-bottom-up"}
#'   or \code{"spectral-top-down"}. \code{"spectral-bottom-up"} grows the
#'   partition one split at a time from a single macro-state, and is the
#'   more reliable choice for a large reduction (\code{k} much smaller than
#'   the number of states); \code{"spectral-top-down"} starts from every
#'   micro-state on its own and repeatedly merges the least costly pair,
#'   which suits a small reduction (\code{k} close to the number of
#'   states). \code{"adaptive"} follows the same rule of thumb as PyDTMC:
#'   top-down below 30 states, otherwise bottom-up when \code{k} is at most
#'   30\% of the number of states and top-down otherwise.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{partition}}{A named list of character vectors giving the
#'       original state names belonging to each macro-state, suitable for
#'       passing to \code{\link{lump}} or \code{\link{is.lumpable}}. Unlike
#'       PyDTMC, which only labels the reduced chain's states generically
#'       (e.g. \code{"ASBU1"}), this traces every macro-state back to the
#'       original states it stands for.}
#'     \item{\code{aggregatedChain}}{The reduced \code{markovchain} object,
#'       row-stochastic, with states named after \code{partition}.}
#'     \item{\code{klDivergence}}{The Kullback-Leibler divergence rate
#'       (in bits) between \code{object} and the lifted \code{aggregatedChain}.
#'       It is zero (up to rounding) when every source state splits its
#'       outgoing probability among a destination macro-state's members in
#'       the same proportions -- those of the stationary distribution
#'       restricted to that macro-state -- regardless of the source; this is
#'       \emph{stronger} than the Kemeny-Snell strong lumpability checked by
#'       \code{\link{is.lumpable}}, which only requires the macro-to-macro
#'       \emph{totals} to agree across sources (see Details).}
#'     \item{\code{method}}{The method actually used, after resolving
#'       \code{"adaptive"}.}
#'     \item{\code{k}}{The number of macro-states actually used, after
#'       resolving an automatic \code{NULL}.}
#'   }
#'
#' @details
#' This targets the same problem as \code{\link{autoLump}}, but by a
#' different and more principled route: \code{autoLump} clusters the
#' leading eigenvectors with k-means (a generic, randomized heuristic,
#' fixed here to a deterministic seed only for reproducibility), while
#' \code{aggregateStates} greedily minimizes the actual information-theoretic
#' quantity that measures how much the aggregation distorts the chain's
#' dynamics.
#'
#' \strong{When is the divergence exactly zero?} Not simply whenever
#' \code{object} is strongly lumpable with respect to \code{partition} in
#' the sense of \code{\link{is.lumpable}}. Strong lumpability only requires
#' that, for every pair of macro-states, all micro-states in the same source
#' macro-state have the same \emph{total} probability of moving to the
#' destination macro-state; it says nothing about how that total is split
#' among the destination macro-state's own members. The Kullback-Leibler
#' divergence used here is sensitive to exactly that split: it is zero (up
#' to rounding) only when every source state distributes its outgoing
#' probability across a destination macro-state's members in the same
#' proportions -- those of the stationary distribution restricted to that
#' macro-state -- regardless of which source state it is. This is a
#' genuinely stronger condition, and a strongly lumpable chain need not
#' satisfy it: the classical Land of Oz weather chain (see the package
#' vignette), lumped into \code{Bad_Weather = \{rainy, snowy\}} and
#' \code{Nice_Weather = \{nice\}}, is strongly lumpable, and
#' \code{aggregateStates} correctly recovers that exact partition as optimal
#' and reproduces \code{\link{lump}}'s aggregated transition matrix, but its
#' divergence is strictly positive, because \code{rainy} and \code{snowy}
#' split their probability between \code{rainy} and \code{snowy} themselves
#' differently from one another.
#'
#' Both methods require the chain to be irreducible (for a unique, strictly
#' positive stationary distribution) and aperiodic (the internal averaging
#' step used by \code{"spectral-top-down"} to re-estimate a working
#' stationary distribution as macro-states are merged assumes convergence,
#' which is not guaranteed for a periodic chain). Use
#' \code{\link{lazyChain}} first to remove periodicity if needed.
#'
#' @references
#' Deng, K., Mehta, P. G. and Meyn, S. P. (2011). Optimal Kullback-Leibler
#' Aggregation via Spectral Theory of Markov Chains. \emph{IEEE Transactions
#' on Automatic Control}, 56(12). \doi{10.1109/TAC.2011.2141350}
#'
#' Fill, J. A. (1991). Eigenvalue bounds on convergence to stationarity for
#' nonreversible Markov chains, with an application to the exclusion
#' process. \emph{The Annals of Applied Probability}, 1(1).
#' \doi{10.1214/aoap/1177005981}
#'
#' @seealso \code{\link{autoLump}}, \code{\link{lump}},
#'   \code{\link{is.lumpable}}, \code{\link{closestReversible}}
#'
#' @examples
#' # A chain aggregated into two macro-states {a,b}/{c,d} where, in addition
#' # to being strongly lumpable, "a" and "b" also split their probability
#' # *within* each destination block identically (0.1/0.1 and 0.4/0.4): this
#' # stronger property is what makes the divergence exactly zero (see
#' # Details for a lumpable-but-nonzero counterexample).
#' statesNames <- c("a", "b", "c", "d")
#' P <- matrix(c(0.1, 0.1, 0.4, 0.4,
#'               0.1, 0.1, 0.4, 0.4,
#'               0.3, 0.3, 0.2, 0.2,
#'               0.3, 0.3, 0.2, 0.2), byrow = TRUE, nrow = 4,
#'             dimnames = list(statesNames, statesNames))
#' mc <- new("markovchain", states = statesNames, transitionMatrix = P)
#'
#' result <- aggregateStates(mc, k = 2)
#' result$partition
#' result$klDivergence # zero up to rounding
#'
#' @exportMethod aggregateStates
setGeneric("aggregateStates", function(object, k = NULL,
                                        method = c("adaptive", "spectral-bottom-up", "spectral-top-down")) {
  standardGeneric("aggregateStates")
})

#' @rdname aggregateStates
setMethod("aggregateStates", "markovchain",
          function(object, k = NULL,
                   method = c("adaptive", "spectral-bottom-up", "spectral-top-down")) {
  method <- match.arg(method)

  if (!is.irreducible(object)) {
    stop("aggregateStates is defined here only for irreducible Markov chains.")
  }
  if (period(object) != 1L) {
    stop(paste0(
      "aggregateStates is defined here only for aperiodic (ergodic) Markov ",
      "chains: the top-down method's working stationary distribution is ",
      "only guaranteed to converge for an aperiodic chain. See lazyChain() ",
      "to remove periodicity first."
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
  if (n < 3L) {
    stop("aggregateStates requires at least 3 states (a 2-state chain cannot be reduced further).")
  }

  stateNames <- states(object)
  pi <- as.numeric(steadyStates(object))
  if (length(pi) != n || any(!is.finite(pi)) || any(pi <= 0)) {
    stop("Unable to obtain a valid, strictly positive stationary distribution.")
  }
  pi <- pi / sum(pi)

  if (is.null(k)) {
    k <- .chooseKByEigengap(P)
  } else {
    if (length(k) != 1L || is.na(k) || !is.numeric(k) || k != as.integer(k)) {
      stop("k must be a single integer.")
    }
    k <- as.integer(k)
    if (k < 2L || k > n - 1L) {
      stop("k must be between 2 and the number of states minus 1.")
    }
  }

  if (method == "adaptive") {
    if (n < 30L) {
      method <- "spectral-top-down"
    } else {
      method <- if ((k / n) <= 0.3) "spectral-bottom-up" else "spectral-top-down"
    }
  }

  result <- if (method == "spectral-bottom-up") {
    .aggregateSpectralBottomUp(P, pi, k)
  } else {
    .aggregateSpectralTopDown(P, pi, k)
  }

  Q <- result$Q
  phi <- result$phi

  if (!isTRUE(all.equal(rowSums(phi), rep(1, n), check.attributes = FALSE))) {
    stop("Internal error: the aggregation partition does not cover every state exactly once.")
  }

  macroStates <- paste0("Macro_", seq_len(k))
  partition <- stats::setNames(vector("list", k), macroStates)
  for (j in seq_len(k)) {
    partition[[j]] <- stateNames[phi[, j] > 0]
  }

  dimnames(Q) <- list(macroStates, macroStates)
  Q <- Q / rowSums(Q)

  klDivergence <- .klDivergenceAggregate(P, pi, phi, Q)

  aggregatedChain <- new("markovchain",
                         states = macroStates,
                         byrow = TRUE,
                         transitionMatrix = Q,
                         name = paste0(object@name, " (aggregated, k=", k, ")"))

  list(
    partition = partition,
    aggregatedChain = aggregatedChain,
    klDivergence = klDivergence,
    method = method,
    k = k
  )
})
