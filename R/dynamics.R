# Internal helper: row-stochastic dense copy of the transition matrix. The
# state order is unchanged by the transpose used for column-stochastic
# storage, so results stay aligned with object@states.
.rowStochasticMatrix <- function(object) {
  P <- as.matrix(object@transitionMatrix)
  if (!object@byrow) {
    P <- t(P)
  }
  if (nrow(P) != ncol(P) || any(!is.finite(P)) || any(P < 0)) {
    stop("The transition matrix must be square, finite, and non-negative.")
  }
  P
}

#' Evolution of a distribution over time
#'
#' Propagates an initial distribution of the states of a discrete-time
#' Markov chain forward in time and returns the whole trajectory of
#' distributions.
#'
#' If \eqn{\mu_0} is the initial distribution (a row vector) and \eqn{P} the
#' row-stochastic transition matrix, the distribution after \eqn{t} steps is
#' \deqn{\mu_t = \mu_{t-1} P = \mu_0 P^t.}
#'
#' @param object A \code{markovchain} object.
#' @param steps A single non-negative whole number: the number of steps to
#'   propagate. With \code{steps = 0} only the initial distribution is
#'   returned.
#' @param initial The initial distribution. Either \code{NULL} (default), for
#'   the uniform distribution over the states; a single state name, for a
#'   point mass on that state; or a numeric vector of non-negative
#'   probabilities summing to one. A named numeric vector is matched to the
#'   states by name, an unnamed one by position.
#' @param lastOnly Logical. If \code{TRUE}, only the distribution after
#'   \code{steps} steps is returned. Defaults to \code{FALSE}.
#'
#' @return If \code{lastOnly = FALSE} (default), a numeric matrix with
#'   \code{steps + 1} rows and one column per state: row \eqn{t} (labelled
#'   \code{"t"}, from \code{"0"}) holds the distribution after \eqn{t}
#'   steps, so the first row is the initial distribution. Otherwise, a named
#'   numeric vector with the distribution after \code{steps} steps.
#'
#' @details
#' The implementation propagates the vector step by step, at a cost of
#' \eqn{O(n^2)} per step for a dense \eqn{n}-state chain, instead of
#' forming \eqn{P^t}; this is what makes the whole trajectory available at
#' no extra cost. Both row- and column-stochastic storage are supported.
#' Each distribution is renormalized after every step to prevent round-off
#' from accumulating over long horizons.
#'
#' This mirrors PyDTMC's \code{redistribute()}, with the same defaults
#' (uniform initial distribution, output including the initial one).
#' For chains that converge, the rows approach the stationary distribution
#' (see \code{\link{steadyStates}}); for periodic chains they do not, which is
#' the expected behaviour and not an error.
#'
#' @seealso \code{\link{steadyStates}}, \code{\link{mixingTime}},
#'   \code{\link{autoplot.markovchain}}
#'
#' @examples
#' statesNames <- c("a", "b")
#' mc <- new("markovchain",
#'   states = statesNames,
#'   transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
#'     byrow = TRUE, nrow = 2,
#'     dimnames = list(statesNames, statesNames)))
#' redistribute(mc, steps = 5, initial = "a")
#' redistribute(mc, steps = 50, initial = c(a = 0.2, b = 0.8), lastOnly = TRUE)
#'
#' @exportMethod redistribute
setGeneric("redistribute", function(object, steps, initial = NULL,
                                    lastOnly = FALSE)
  standardGeneric("redistribute"))

#' @rdname redistribute
setMethod("redistribute", "markovchain",
  function(object, steps, initial = NULL, lastOnly = FALSE) {
    if (length(steps) != 1L || !is.numeric(steps) || !is.finite(steps) ||
        steps < 0 || steps != round(steps)) {
      stop("steps must be a single non-negative whole number.")
    }
    if (length(lastOnly) != 1L || !is.logical(lastOnly) || is.na(lastOnly)) {
      stop("lastOnly must be TRUE or FALSE.")
    }
    steps <- as.integer(steps)

    P <- .rowStochasticMatrix(object)
    states <- object@states
    n <- length(states)

    if (is.null(initial)) {
      mu <- rep(1 / n, n)
    } else if (is.character(initial)) {
      if (length(initial) != 1L || !(initial %in% states)) {
        stop("A character initial must be a single state of the chain.")
      }
      mu <- as.numeric(states == initial)
    } else if (is.numeric(initial)) {
      if (length(initial) != n || any(!is.finite(initial)) ||
          any(initial < 0)) {
        stop(paste0(
          "A numeric initial must contain ", n,
          " finite, non-negative probabilities."
        ))
      }
      if (!is.null(names(initial))) {
        if (anyDuplicated(names(initial)) ||
            !setequal(names(initial), states)) {
          stop("Names of initial must match the states of the chain.")
        }
        initial <- initial[states]
      }
      if (abs(sum(initial) - 1) > sqrt(.Machine$double.eps) * n) {
        stop("initial must sum to one.")
      }
      mu <- as.numeric(initial) / sum(initial)
    } else {
      stop("initial must be NULL, a state name, or a numeric vector.")
    }

    trajectory <- matrix(0, nrow = steps + 1L, ncol = n,
                         dimnames = list(as.character(0:steps), states))
    trajectory[1L, ] <- mu
    for (t in seq_len(steps)) {
      nextMu <- as.numeric(trajectory[t, ] %*% P)
      trajectory[t + 1L, ] <- nextMu / sum(nextMu)
    }

    if (lastOnly) {
      return(trajectory[steps + 1L, ])
    }
    trajectory
  })

# Internal helper: spectral radius (Perron root) of the 0/1 adjacency matrix
# of the transition graph. Always >= 1, since every row of a stochastic
# matrix has a positive entry, so every state leads to a cycle.
.adjacencyPerronRoot <- function(object) {
  A <- (.rowStochasticMatrix(object) > 0) * 1
  max(Mod(eigen(A, only.values = TRUE)$values))
}

#' Topological entropy of a Markov chain
#'
#' Computes the topological entropy of the graph of a discrete-time Markov
#' chain: the exponential growth rate of the number of distinct admissible
#' paths, ignoring their probabilities.
#'
#' If \eqn{A} is the 0/1 adjacency matrix with \eqn{A_{ij}=1} exactly when
#' \eqn{p_{ij}>0}, the topological entropy is
#' \deqn{h_{top} = \log_b \rho(A),}
#' where \eqn{\rho(A)} is the spectral radius (Perron root) of \eqn{A}.
#'
#' @param object A \code{markovchain} object.
#' @param base A finite numeric scalar strictly greater than one. The default,
#'   \code{2}, returns bits, consistently with \code{\link{entropyRate}}.
#'
#' @return A non-negative numeric scalar in units determined by \code{base}.
#'
#' @details
#' Only the pattern of positive entries matters: the value depends on which
#' transitions are possible, not on how likely they are. It is the upper
#' bound of the entropy rate over all the Markov chains sharing that graph
#' (the variational principle, see Parry, 1964), so
#' \code{entropyRate(object) <= topologicalEntropy(object)} for an
#' irreducible chain. The bound is attained by the maximal-entropy
#' (Parry) chain on the same graph, and also, for instance, by a chain whose
#' every row is uniform over a common number of successors. A chain that is a
#' single cycle (deterministic dynamics) has \code{topologicalEntropy = 0}.
#'
#' No irreducibility is needed: for a reducible chain the result is the
#' largest value over its communicating classes. A probability that is
#' positive but numerically tiny counts as a transition, exactly as in
#' \code{\link{is.irreducible}}.
#'
#' The cost is one eigenvalue computation, \eqn{O(n^3)} time and
#' \eqn{O(n^2)} memory for a dense chain. It mirrors PyDTMC's
#' \code{topological_entropy}, which uses the natural logarithm.
#'
#' @references
#' Parry, W. (1964). Intrinsic Markov chains. \emph{Transactions of the
#' American Mathematical Society}, 112, 55-66.
#'
#' Cover, T. M. and Thomas, J. A. (2006). \emph{Elements of Information
#' Theory}, 2nd edition. Wiley.
#'
#' @seealso \code{\link{entropyRate}}, \code{\link{normalizedEntropyRate}}
#'
#' @examples
#' statesNames <- c("a", "b")
#' mc <- new("markovchain",
#'   states = statesNames,
#'   transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
#'     byrow = TRUE, nrow = 2,
#'     dimnames = list(statesNames, statesNames)))
#' topologicalEntropy(mc)
#'
#' @exportMethod topologicalEntropy
setGeneric("topologicalEntropy", function(object, base = 2)
  standardGeneric("topologicalEntropy"))

#' @rdname topologicalEntropy
setMethod("topologicalEntropy", "markovchain", function(object, base = 2) {
  if (length(base) != 1L || !is.numeric(base) || !is.finite(base) ||
      base <= 1) {
    stop("base must be a single finite number strictly greater than 1.")
  }
  value <- log(.adjacencyPerronRoot(object)) / log(base)
  if (!is.finite(value)) {
    stop("Unable to compute a finite topological entropy.")
  }
  # The Perron root is at least one; discard round-off around zero.
  if (value < 100 * .Machine$double.eps) {
    value <- 0
  }
  as.numeric(value)
})

#' Normalized entropy rate of a Markov chain
#'
#' The entropy rate of a finite irreducible Markov chain divided by the
#' topological entropy of its graph, a measure of how random the chain is
#' relative to the most random chain with the same possible transitions.
#'
#' @param object A \code{markovchain} object representing a finite,
#'   irreducible discrete-time Markov chain.
#'
#' @return A numeric scalar in \eqn{[0,1]}: \code{0} for deterministic
#'   dynamics, \code{1} when the entropy rate reaches the topological
#'   entropy.
#'
#' @details
#' The ratio \eqn{H / h_{top}} does not depend on the logarithm base, so no
#' \code{base} argument is needed. By the variational principle
#' (\code{\link{topologicalEntropy}}) it never exceeds one; tiny excursions
#' due to round-off are clipped to \eqn{[0,1]}. When \eqn{h_{top}=0} the
#' chain has a single possible path, hence entropy rate zero, and the
#' ratio is defined here to be \code{0}, the convention used by PyDTMC's
#' \code{entropy_rate_normalized}.
#'
#' Irreducibility is required, as in \code{\link{entropyRate}}, which is
#' what makes the entropy rate well defined.
#'
#' @seealso \code{\link{entropyRate}}, \code{\link{topologicalEntropy}}
#'
#' @examples
#' statesNames <- c("a", "b")
#' mc <- new("markovchain",
#'   states = statesNames,
#'   transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
#'     byrow = TRUE, nrow = 2,
#'     dimnames = list(statesNames, statesNames)))
#' normalizedEntropyRate(mc)
#'
#' @exportMethod normalizedEntropyRate
setGeneric("normalizedEntropyRate", function(object)
  standardGeneric("normalizedEntropyRate"))

#' @rdname normalizedEntropyRate
setMethod("normalizedEntropyRate", "markovchain", function(object) {
  h <- entropyRate(object, base = 2)
  htop <- topologicalEntropy(object, base = 2)
  if (h <= 100 * .Machine$double.eps || htop <= 100 * .Machine$double.eps) {
    return(0)
  }
  as.numeric(min(1, max(0, h / htop)))
})

#' Relaxation time of a Markov chain
#'
#' The relaxation time \eqn{1/\mathrm{gap}} of a finite irreducible
#' discrete-time Markov chain, the reciprocal of its spectral gap.
#'
#' @param object A \code{markovchain} object representing a finite,
#'   irreducible discrete-time Markov chain.
#'
#' @return A positive numeric scalar, in number of steps: \code{Inf} for a
#'   periodic chain (spectral gap zero), \code{1} for the trivial one-state
#'   chain.
#'
#' @details
#' Following Levin and Peres (2017, Section 12.2) the relaxation time is
#' \eqn{t_{rel} = 1/\gamma}, with \eqn{\gamma = 1 - \mathrm{SLEM}} the
#' spectral gap of \code{\link{spectralGap}}. It is the quantity PyDTMC
#' calls \code{relaxation_rate}, even if it is a time, not a rate. The
#' two are the same number.
#'
#' PyDTMC's \code{mixing_rate} is \eqn{-1/\log(\mathrm{SLEM})}, i.e. the
#' timescale of the slowest non-trivial mode, which is already returned
#' as the first element (\code{"tau2"}) of \code{\link{impliedTimescales}}.
#' It is therefore not repeated as a separate function.
#'
#' @references
#' Levin, D. A. and Peres, Y. (2017). \emph{Markov Chains and Mixing Times},
#' 2nd edition. American Mathematical Society.
#'
#' @seealso \code{\link{spectralGap}}, \code{\link{slem}},
#'   \code{\link{impliedTimescales}}, \code{\link{mixingTime}}
#'
#' @examples
#' statesNames <- c("a", "b")
#' mc <- new("markovchain",
#'   states = statesNames,
#'   transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
#'     byrow = TRUE, nrow = 2,
#'     dimnames = list(statesNames, statesNames)))
#' relaxationTime(mc)
#'
#' @exportMethod relaxationTime
setGeneric("relaxationTime", function(object)
  standardGeneric("relaxationTime"))

#' @rdname relaxationTime
setMethod("relaxationTime", "markovchain", function(object) {
  gap <- spectralGap(object)
  if (gap <= sqrt(.Machine$double.eps)) {
    return(Inf)
  }
  as.numeric(1 / gap)
})
