# Distances, divergences and a test of equality between two Markov chains on
# the same states.
#
# For a junior reader: two transition matrices P and Q can differ in many
# ways, and no single number says everything. Row-wise measures (total
# variation, Hellinger) compare the next-step distribution from each state;
# the Kullback-Leibler divergence *rate* weights those row differences by how
# often the chain visits each state, which is what matters for long runs of
# the process; the distance between stationary distributions only looks at
# the long run. When the chains are estimated from data, a test tells
# whether the differences are larger than sampling noise.

# Internal helper: the transition matrix of a chain as a row-stochastic
# matrix with its states as dimnames (byrow = FALSE chains are transposed).
.rowStochastic <- function(object) {
  P <- as.matrix(object@transitionMatrix)
  if (!object@byrow) {
    P <- t(P)
  }
  dimnames(P) <- list(object@states, object@states)
  P
}

# Internal helper: Kullback-Leibler divergence of each row of P from the same
# row of Q, with 0 log 0 = 0 and Inf where P > 0 and Q = 0.
.rowKL <- function(P, Q, base = exp(1)) {
  terms <- P * log(P / Q)
  terms[P == 0] <- 0
  rowSums(terms) / log(base)
}

# Internal helper: check that two chains share their states and return their
# row-stochastic matrices, the second in the order of the first.
.alignChains <- function(x, y) {
  if (!is(x, "markovchain") || !is(y, "markovchain")) {
    stop("x and y must be markovchain objects.")
  }
  if (!setequal(x@states, y@states) || length(x@states) != length(y@states)) {
    stop("x and y must have the same states (they are matched by name).")
  }
  P <- .rowStochastic(x)
  Q <- .rowStochastic(y)[x@states, x@states, drop = FALSE]
  list(P = P, Q = Q)
}

# Internal helper: the stationary distributions of a chain, one per row
# (steadyStates() returns them as columns for column-stochastic chains).
.stationaryRows <- function(object) {
  ss <- steadyStates(object)
  if (!object@byrow) {
    ss <- t(ss)
  }
  ss
}

# Internal helper: the unique stationary distribution of a chain, or an
# informative error.
.uniqueStationary <- function(object, argName) {
  ss <- .stationaryRows(object)
  if (nrow(ss) != 1L) {
    stop(argName, " has ", nrow(ss), " stationary distributions; ",
         "a unique one is needed (supply weights to use other ones).")
  }
  stats::setNames(as.numeric(ss), object@states)
}

#' Kullback-Leibler divergence rate between two Markov chains
#'
#' Computes the Kullback-Leibler divergence rate of a Markov chain \code{x}
#' from a Markov chain \code{y} on the same states: the limit of
#' \eqn{\frac{1}{n} D(X_1,\ldots,X_n \,\|\, Y_1,\ldots,Y_n)}, the average
#' information per step for discriminating \code{x} against \code{y}.
#'
#' @param x,y Two \code{markovchain} objects with the same states, matched
#'   by name.
#' @param weights Optional non-negative numeric vector, one entry per state
#'   (named or in the order of \code{states(x)}), with which the row
#'   divergences are averaged; it is normalized to sum to one. The default,
#'   \code{NULL}, uses the stationary distribution of \code{x}, which must
#'   then be unique.
#' @param base The base of the logarithm: the default \code{exp(1)} gives the
#'   rate in nats, \code{2} in bits.
#'
#' @return A single non-negative number, \code{Inf} when \code{x} makes with
#'   positive probability a transition that \code{y} forbids from a state of
#'   positive weight.
#'
#' @details
#' For chains with transition matrices \eqn{P} (of \code{x}) and \eqn{Q}
#' (of \code{y}) and with \eqn{\pi} the stationary distribution of \eqn{P},
#' the divergence rate is
#' \deqn{D(P \,\|\, Q) = \sum_i \pi_i \sum_j P_{ij} \log \frac{P_{ij}}{Q_{ij}},}
#' which Rached, Alajaji and Campbell (2004) prove to be the limit above
#' for irreducible chains, whatever the initial distributions. By the
#' ergodic theorem it is also the almost sure limit of
#' \eqn{\frac1n \log \frac{\Pr_P(X_1, \ldots, X_n)}{\Pr_Q(X_1, \ldots, X_n)}}
#' along a path of \code{x}. It is zero exactly when the rows of \eqn{P} and
#' \eqn{Q} agree on every state of positive weight, and it is not symmetric.
#'
#' With \code{weights} equal to the relative frequencies of the departure
#' states of a sequence, \eqn{2N D(\hat P \,\|\, Q)}, with \eqn{N} the number
#' of transitions, is the likelihood-ratio statistic \eqn{G^2} of the
#' hypothesis that the sequence comes from \eqn{Q} (Kullback, Kupperman and
#' Ku 1962), as computed by \code{\link{verifyEmpiricalToTheoretical}}.
#'
#' @references
#' Rached, Z., Alajaji, F. and Campbell, L. L. (2004). The Kullback-Leibler
#' divergence rate between Markov sources. \emph{IEEE Transactions on
#' Information Theory}, 50(5), 917-921. \doi{10.1109/TIT.2004.826687}
#'
#' Kullback, S., Kupperman, M. and Ku, H. H. (1962). Tests for contingency
#' tables and Markov chains. \emph{Technometrics}, 4(4), 573-608.
#'
#' @seealso \code{\link{compareChains}}, \code{\link{entropyRate}}
#'
#' @examples
#' s <- c("a", "b", "c")
#' P <- new("markovchain", states = s, transitionMatrix = matrix(
#'   c(0.5, 0.3, 0.2, 0.1, 0.6, 0.3, 0.25, 0.25, 0.5), 3, byrow = TRUE,
#'   dimnames = list(s, s)))
#' Q <- new("markovchain", states = s, transitionMatrix = matrix(
#'   c(0.4, 0.4, 0.2, 0.2, 0.5, 0.3, 0.3, 0.1, 0.6), 3, byrow = TRUE,
#'   dimnames = list(s, s)))
#' divergenceRate(P, Q)
#' divergenceRate(Q, P) # not symmetric
#' divergenceRate(P, P)
#'
#' @export
divergenceRate <- function(x, y, weights = NULL, base = exp(1)) {
  if (length(base) != 1L || !is.numeric(base) || !is.finite(base) || base <= 0 || base == 1) {
    stop("base must be a single finite positive number other than 1.")
  }
  m <- .alignChains(x, y)
  if (is.null(weights)) {
    w <- .uniqueStationary(x, "x")
  } else {
    w <- .stateWeights(weights, x@states)
  }
  rows <- .rowKL(m$P, m$Q, base)
  # states of zero weight do not contribute, even with an infinite row
  sum(w[w > 0] * rows[w > 0])
}

# Internal helper: validate a weight vector over the states.
.stateWeights <- function(weights, stateNames) {
  if (!is.numeric(weights) || length(weights) != length(stateNames) ||
      any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
    stop("weights must be a finite non-negative vector with one entry per state ",
         "and a positive sum.")
  }
  if (!is.null(names(weights))) {
    if (!setequal(names(weights), stateNames)) {
      stop("the names of weights must be the states.")
    }
    weights <- weights[stateNames]
  }
  stats::setNames(as.numeric(weights) / sum(weights), stateNames)
}

#' Compare two Markov chains
#'
#' Computes distances and divergences between two Markov chains on the same
#' states and, when the chains are given as observed sequences, tests
#' whether the sequences come from the same chain.
#'
#' @param x,y Two \code{markovchain} objects with the same states (matched by
#'   name), or two observed sequences (character, factor or numeric
#'   vectors), from which the chains are estimated by maximum likelihood on
#'   the union of their states.
#' @param base The base of the logarithm used by the divergences (default
#'   \code{exp(1)}, nats).
#' @param test For sequences, the test of equality of the two transition
#'   matrices passed to \code{\link{verifyHomogeneity}}: \code{"G"}
#'   (likelihood ratio, the default), \code{"Pearson"}, \code{"simulation"}
#'   (Monte Carlo), or \code{"none"}. Ignored for \code{markovchain}
#'   objects, whose sampling variability is not known.
#' @param ... Further arguments for \code{\link{verifyHomogeneity}}, such as
#'   \code{B} and \code{seed}.
#'
#' @return An object of class \code{"markovchainComparison"}, a list with:
#'   \describe{
#'     \item{\code{distances}}{A named vector: \code{totalVariation}, the
#'       largest row-wise total variation distance
#'       \eqn{\max_i \frac12\sum_j |P_{ij} - Q_{ij}|};
#'       \code{weightedTotalVariation}, the same row distances averaged with
#'       the stationary distribution of \code{x}; \code{hellinger}, the
#'       largest row-wise Hellinger distance
#'       \eqn{\max_i \sqrt{\frac12 \sum_j (\sqrt{P_{ij}} - \sqrt{Q_{ij}})^2}};
#'       \code{frobenius}, \eqn{\|P - Q\|_F}; \code{maximum},
#'       \eqn{\max_{i,j}|P_{ij} - Q_{ij}|}; and \code{stationary}, the total
#'       variation distance between the two stationary distributions.}
#'     \item{\code{divergences}}{A named vector: \code{rateXY} and
#'       \code{rateYX}, the Kullback-Leibler divergence rates
#'       \eqn{D(P\|Q)} and \eqn{D(Q\|P)} (see \code{\link{divergenceRate}});
#'       \code{maxRowXY}, the largest row-wise Kullback-Leibler divergence
#'       of \eqn{P} from \eqn{Q}; and \code{jensenShannon}, the largest
#'       row-wise Jensen-Shannon divergence, which is symmetric and finite.}
#'     \item{\code{stationaryDistributions}}{A two-row matrix with the
#'       stationary distributions of \code{x} and \code{y}.}
#'     \item{\code{test}}{For sequences, the \code{htest} returned by
#'       \code{\link{verifyHomogeneity}}; otherwise \code{NULL}.}
#'     \item{\code{x}, \code{y}}{The two chains, as \code{markovchain}
#'       objects.}
#'   }
#'
#' @details
#' The measures that depend on a stationary distribution
#' (\code{weightedTotalVariation}, \code{stationary}, \code{rateXY},
#' \code{rateYX}) are \code{NA} when the chain concerned has more than one.
#' The row-wise measures and their names follow PyDTMC's
#' \code{calculate_transition_distance},
#' \code{calculate_transition_divergence} and
#' \code{calculate_stationary_distance}, which give the same values.
#'
#' For two sequences the test is that of Anderson and Goodman (1957) and
#' Kullback, Kupperman and Ku (1962): the hypothesis that both come from the
#' same chain, against different transition matrices. Each sequence is
#' treated as one sample; for several sequences per group, or for counts,
#' call \code{\link{verifyHomogeneity}} directly.
#'
#' @references
#' Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about
#' Markov chains. \emph{The Annals of Mathematical Statistics}, 28(1),
#' 89-110.
#'
#' Kullback, S., Kupperman, M. and Ku, H. H. (1962). Tests for contingency
#' tables and Markov chains. \emph{Technometrics}, 4(4), 573-608.
#'
#' Rached, Z., Alajaji, F. and Campbell, L. L. (2004). The Kullback-Leibler
#' divergence rate between Markov sources. \emph{IEEE Transactions on
#' Information Theory}, 50(5), 917-921. \doi{10.1109/TIT.2004.826687}
#'
#' @seealso \code{\link{divergenceRate}}, \code{\link{verifyHomogeneity}},
#'   \code{\link{autoplot.markovchain}} (\code{type = "comparison"})
#'
#' @examples
#' s <- c("a", "b", "c")
#' P <- new("markovchain", states = s, transitionMatrix = matrix(
#'   c(0.5, 0.3, 0.2, 0.1, 0.6, 0.3, 0.25, 0.25, 0.5), 3, byrow = TRUE,
#'   dimnames = list(s, s)))
#' Q <- new("markovchain", states = s, transitionMatrix = matrix(
#'   c(0.4, 0.4, 0.2, 0.2, 0.5, 0.3, 0.3, 0.1, 0.6), 3, byrow = TRUE,
#'   dimnames = list(s, s)))
#' compareChains(P, Q)
#'
#' # two sequences: distances between the estimated chains, and a test
#' set.seed(1)
#' x <- rmarkovchain(300, P, t0 = "a")
#' y <- rmarkovchain(300, Q, t0 = "a")
#' compareChains(x, y)
#'
#' @export
compareChains <- function(x, y, base = exp(1), test = c("G", "Pearson", "simulation", "none"),
                          ...) {
  test <- match.arg(test)
  testResult <- NULL
  isChain <- function(z) is(z, "markovchain")
  if (!isChain(x) || !isChain(y)) {
    # Step 1 (sequences): estimate both chains on the same states.
    isSequence <- function(z) is.atomic(z) && !is.matrix(z) && length(z) >= 2L
    if (!isSequence(x) || !isSequence(y)) {
      stop("x and y must both be markovchain objects or both be sequences.")
    }
    xs <- as.character(x)
    ys <- as.character(y)
    if (anyNA(xs) || anyNA(ys)) {
      stop("the sequences must not contain missing values.")
    }
    allStates <- sort(unique(c(xs, ys)))
    # estimate(): maximum likelihood chain (row-normalized counts) on the common set of states
    estimate <- function(seq, label) {
      counts <- createSequenceMatrix(seq, possibleStates = allStates)
      counts <- counts[allStates, allStates, drop = FALSE]
      never <- allStates[rowSums(counts) == 0]
      if (length(never)) {
        stop("in ", label, " the state(s) ", paste(never, collapse = ", "),
             " are never left, so their transitions cannot be estimated; ",
             "fit the chains yourself and pass markovchain objects.")
      }
      new("markovchain", states = allStates, byrow = TRUE,
          transitionMatrix = counts / rowSums(counts), name = label)
    }
    xChain <- estimate(xs, "x")
    yChain <- estimate(ys, "y")
    # optional test that the two sequences share the same transition matrix
    if (test != "none") {
      testResult <- verifyHomogeneity(list(x = xs, y = ys), method = test,
                                      verbose = FALSE, ...)
    }
  } else {
    xChain <- x
    yChain <- y
  }

  if (length(base) != 1L || !is.numeric(base) || !is.finite(base) || base <= 0 || base == 1) {
    stop("base must be a single finite positive number other than 1.")
  }
  # Step 2: the matrices, in the order of the states of x.
  m <- .alignChains(xChain, yChain)
  P <- m$P
  Q <- m$Q
  # stationaryOrNA(): the stationary distribution in the order of the states of x, or NA
  # for every state if it is not unique (the weighted summaries are then not defined)
  stationaryOrNA <- function(chain) {
    ss <- .stationaryRows(chain)
    if (nrow(ss) == 1L) {
      stats::setNames(as.numeric(ss), chain@states)[xChain@states]
    } else {
      stats::setNames(rep(NA_real_, length(xChain@states)), xChain@states)
    }
  }
  piX <- stationaryOrNA(xChain)
  piY <- stationaryOrNA(yChain)

  # Step 3: row-wise distances and their summaries. Row i of P and of Q are the two
  # distributions of the next state given the current state i.
  rowTV <- 0.5 * rowSums(abs(P - Q))
  rowHellinger <- sqrt(0.5 * rowSums((sqrt(P) - sqrt(Q))^2))
  distances <- c(
    totalVariation = max(rowTV),            # worst row
    # (the other entries: average over the stationary law of x, worst Hellinger row,
    # Frobenius norm and largest entry of P - Q, total variation of the stationary laws)
    weightedTotalVariation = sum(piX * rowTV),
    hellinger = max(rowHellinger),
    frobenius = sqrt(sum((P - Q)^2)),
    maximum = max(abs(P - Q)),
    stationary = 0.5 * sum(abs(piX - piY))
  )

  # Step 4: divergences. The rates weight the row divergences by the
  # stationary distribution of the first chain; states of weight zero
  # contribute nothing even if their row divergence is infinite.
  rowXY <- .rowKL(P, Q, base)        # KL divergence of row i of P from row i of Q
  rowYX <- .rowKL(Q, P, base)
  M <- (P + Q) / 2                   # mixture used by the Jensen-Shannon divergence
  rowJS <- 0.5 * (.rowKL(P, M, base) + .rowKL(Q, M, base))
  # rate(): weighted average of the row divergences, NA if the weights are not available
  rate <- function(w, rows) if (anyNA(w)) NA_real_ else sum(w[w > 0] * rows[w > 0])
  divergences <- c(
    rateXY = rate(piX, rowXY),
    rateYX = rate(piY, rowYX),
    maxRowXY = max(rowXY),
    jensenShannon = max(rowJS)
  )

  structure(list(distances = distances, divergences = divergences,
                 stationaryDistributions = rbind(x = piX, y = piY),
                 test = testResult, x = xChain, y = yChain),
            class = "markovchainComparison")
}

#' Print a comparison of two Markov chains
#'
#' @param x A \code{"markovchainComparison"} object returned by
#'   \code{\link{compareChains}}.
#' @param digits Number of significant digits.
#' @param ... Passed to \code{print}.
#' @return \code{x}, invisibly.
#' @method print markovchainComparison
#' @keywords internal
#' @export
print.markovchainComparison <- function(x, digits = 4, ...) {
  cat("Comparison of two Markov chains on", length(x$x@states), "states\n\n")
  cat("Distances:\n")
  print(signif(x$distances, digits), ...)
  cat("\nDivergences (rates weighted by the stationary distributions):\n")
  print(signif(x$divergences, digits), ...)
  if (!is.null(x$test)) {
    cat("\nTest that both sequences come from the same chain:\n")
    cat(sprintf("  %s = %s, df = %s, p-value = %s\n", names(x$test$statistic),
                format(signif(x$test$statistic, digits)), format(x$test$parameter),
                format.pval(x$test$p.value, digits = digits)))
  }
  invisible(x)
}
