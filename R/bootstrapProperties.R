# Bootstrap confidence intervals for properties of a Markov chain estimated
# from data (SLEM, spectral gap, relaxation time, mixing time, ...).
#
# For a junior reader: the functions slem(), spectralGap(), mixingTime() and
# the like return a single number for a chain. When the chain is estimated
# from a sequence, that number inherits the sampling error of the estimated
# transition matrix, and the error is not easy to derive by hand (these are
# non-linear functions of the matrix, some not even smooth). The bootstrap
# approximates it by re-estimating the chain on many artificial samples drawn
# from the fitted one, computing the property each time, and reading the
# interval off the spread of the results.

# Internal helper: the built-in statistics, as functions of a markovchain.
.bootstrapBuiltins <- function(epsilon) {
  list(
    slem = function(mc) slem(mc),
    spectralGap = function(mc) spectralGap(mc),
    relaxationTime = function(mc) relaxationTime(mc),
    mixingTime = function(mc) mixingTime(mc, epsilon = epsilon),
    entropyRate = function(mc) entropyRate(mc),
    kemenyConstant = function(mc) kemenyConstant(mc),
    steadyStates = function(mc) {
      ss <- steadyStates(mc)
      if (nrow(ss) != 1L) stop("more than one stationary distribution")
      stats::setNames(as.numeric(ss), mc@states)
    }
  )
}

# Internal helper: evaluate every statistic on one chain and flatten the
# results into a named vector. A statistic that fails (for example a
# replicate that is not irreducible) gives NA values of the right length.
.evaluateStatistics <- function(chain, statistics, widths) {
  unlist(lapply(names(statistics), function(nm) {
    value <- tryCatch(as.numeric(statistics[[nm]](chain)), error = function(e) NULL)
    if (is.null(value) || length(value) != widths[[nm]]) {
      value <- rep(NA_real_, widths[[nm]])
    }
    value
  }), use.names = FALSE)
}

#' Bootstrap confidence intervals for properties of an estimated chain
#'
#' Estimates a Markov chain from one or more observed sequences and gives
#' bootstrap standard errors and confidence intervals for properties of the
#' chain, such as the second largest eigenvalue modulus, the spectral gap,
#' the relaxation time or the mixing time, which the corresponding functions
#' return as point values only.
#'
#' @param data An observed sequence (character, factor or numeric vector),
#'   a list of independent sequences, or a square matrix of transition counts
#'   with the states as row and column names (only with
#'   \code{method = "multinomial"}).
#' @param statistics Either a character vector of built-in statistics, among
#'   \code{"slem"}, \code{"spectralGap"}, \code{"relaxationTime"},
#'   \code{"mixingTime"} (the default four), \code{"entropyRate"},
#'   \code{"kemenyConstant"} and \code{"steadyStates"}, or a named list
#'   whose elements are such names or functions of a \code{markovchain}
#'   object returning a numeric vector of fixed length (an unnamed built-in
#'   name is named after itself).
#' @param nboot A single integer of at least \code{2}: the number of
#'   bootstrap replicates.
#' @param method How the replicates are drawn: \code{"multinomial"} (the
#'   default) resamples, for each state, its observed number of departures
#'   from the estimated row; \code{"parametric"} simulates new sequences of
#'   the observed lengths from the estimated chain, each starting from the
#'   observed first state. See Details.
#' @param level The confidence level, a single number in \eqn{(0, 1)}.
#' @param type The kind of interval: \code{"percentile"} (the default) or
#'   \code{"basic"}.
#' @param possibleStates Optional character vector of states to include
#'   even if they are not observed. Every state must be left at least once in
#'   the data, otherwise its transition probabilities cannot be estimated.
#' @param epsilon The total variation threshold of \code{\link{mixingTime}}.
#' @param seed An optional whole number. When given, the replicates are drawn
#'   after \code{set.seed(seed)} and the caller's random number stream is left
#'   untouched.
#'
#' @return An object of class \code{"markovchainBootstrap"}, a list with:
#'   \describe{
#'     \item{\code{estimate}}{The statistics of the estimated chain.}
#'     \item{\code{ci}}{A matrix with one row per statistic and the lower and
#'       upper limits of the interval.}
#'     \item{\code{se}, \code{bias}}{Bootstrap standard errors and bias
#'       estimates (mean of the replicates minus the estimate).}
#'     \item{\code{failures}}{For each statistic, the number of replicates
#'       for which it could not be computed (they are left out of the
#'       intervals).}
#'     \item{\code{replicates}}{The \code{nboot} by statistics matrix of
#'       bootstrap values.}
#'     \item{\code{chain}}{The estimated \code{markovchain}.}
#'     \item{\code{counts}}{The matrix of observed transition counts.}
#'     \item{\code{method}, \code{type}, \code{level}, \code{nboot}}{The
#'       settings used.}
#'   }
#'
#' @details
#' The chain is estimated by maximum likelihood, \eqn{\hat P_{ij} =
#' n_{ij}/n_i}. With \code{method = "multinomial"} a replicate draws, for
#' every state \eqn{i}, the counts \eqn{(n^*_{i1}, \ldots, n^*_{ir})} from a
#' multinomial distribution with \eqn{n_i} trials and probabilities
#' \eqn{\hat P_{i\cdot}}, and re-estimates \eqn{P^*_{ij} = n^*_{ij}/n_i}. This
#' is the conditional bootstrap of Basawa et al. (1990), valid for finite
#' ergodic chains, which keeps the number of visits to each state fixed. With
#' \code{method = "parametric"} a replicate simulates a new path of each
#' observed length from \eqn{\hat P}, as in Kulperger and Prakasa Rao
#' (1989); the numbers of visits then vary too, and a replicate in which some
#' state is never left cannot be estimated and counts as a failure. The
#' parametric version is slower but closer to the actual sampling scheme of a
#' single path.
#'
#' The percentile interval is formed by the \eqn{(1-\code{level})/2} and
#' \eqn{(1+\code{level})/2} quantiles of the replicates; the basic interval
#' reflects them around the estimate, \eqn{2\hat\theta - q_{1-\alpha/2}} and
#' \eqn{2\hat\theta - q_{\alpha/2}}.
#'
#' Some care is needed in reading the results. The SLEM and the spectral gap
#' are not smooth functions of the matrix where the second largest
#' eigenvalue modulus is reached by more than one eigenvalue, the mixing time
#' is an integer, and replicates of a chain with rare transitions may be
#' reducible or periodic, in which case some statistics fail; the bootstrap
#' distribution of all of them can be skewed. The intervals are therefore
#' approximate, and their coverage should be checked by simulation for the
#' case at hand, as done in the package's vignette. The asymptotic theory of
#' these estimators is developed by Hsu et al. (2019) for the spectral gap of
#' reversible chains.
#'
#' @references
#' Basawa, I. V., Green, T. A., McCormick, W. P. and Taylor, R. L. (1990).
#' Asymptotic bootstrap validity for finite Markov chains.
#' \emph{Communications in Statistics - Theory and Methods}, 19(4),
#' 1493-1510.
#'
#' Kulperger, R. J. and Prakasa Rao, B. L. S. (1989). Bootstrapping a finite
#' state Markov chain. \emph{Sankhya, Series A}, 51(2), 178-191.
#'
#' Hsu, D., Kontorovich, A., Levin, D. A., Peres, Y., Szepesvari, C. and
#' Wolfer, G. (2019). Mixing time estimation in reversible Markov chains from
#' a single sample path. \emph{The Annals of Applied Probability}, 29(4),
#' 2439-2480. \doi{10.1214/18-AAP1457}
#'
#' @seealso \code{\link{slem}}, \code{\link{spectralGap}},
#'   \code{\link{relaxationTime}}, \code{\link{mixingTime}},
#'   \code{\link{markovchainFit}} (whose \code{method = "bootstrap"} gives
#'   bootstrap standard errors of the transition probabilities themselves)
#'
#' @examples
#' data(rain)
#' set.seed(1)
#' boot <- bootstrapProperties(rain$rain, nboot = 200)
#' boot
#' boot$ci
#'
#' # any function of the chain, here the stationary probability of no rain
#' bootstrapProperties(rain$rain, nboot = 200, seed = 2,
#'   statistics = list(dry = function(mc) steadyStates(mc)[1, "0"]))
#'
#' @export
bootstrapProperties <- function(data, statistics = c("slem", "spectralGap", "relaxationTime",
                                                     "mixingTime"),
                                nboot = 1000L, method = c("multinomial", "parametric"),
                                level = 0.95, type = c("percentile", "basic"),
                                possibleStates = character(), epsilon = 0.25, seed = NULL) {
  # Step 1: check the arguments.
  method <- match.arg(method)
  type <- match.arg(type)
  if (length(nboot) != 1L || !is.numeric(nboot) || is.na(nboot) ||
      nboot != as.integer(nboot) || nboot < 2L) {
    stop("nboot must be a single integer of at least 2.")
  }
  nboot <- as.integer(nboot)
  if (length(level) != 1L || !is.numeric(level) || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be a single number strictly between 0 and 1.")
  }
  builtins <- .bootstrapBuiltins(epsilon)
  if (is.character(statistics)) {
    unknown <- setdiff(statistics, names(builtins))
    if (length(statistics) == 0L || length(unknown)) {
      stop("unknown statistics: ", paste(unknown, collapse = ", "),
           ". Use some of ", paste(names(builtins), collapse = ", "),
           " or a named list of functions.")
    }
    statistics <- builtins[unique(statistics)]
  } else if (is.list(statistics)) {
    # an unnamed built-in, as in list(dry = f, "slem"), is named after itself
    if (length(statistics)) {
      nms <- names(statistics)
      if (is.null(nms)) nms <- rep("", length(statistics))
      byName <- !nzchar(nms) & vapply(statistics, function(s) is.character(s) && length(s) == 1L,
                                      logical(1))
      nms[byName] <- unlist(statistics[byName])
      names(statistics) <- nms
    }
    if (length(statistics) == 0L || is.null(names(statistics)) ||
        any(!nzchar(names(statistics))) || anyDuplicated(names(statistics))) {
      stop("a list of statistics must have distinct, non-empty names.")
    }
    statistics <- lapply(statistics, function(s) {
      if (is.function(s)) {
        s
      } else if (is.character(s) && length(s) == 1L && s %in% names(builtins)) {
        builtins[[s]]
      } else {
        stop("each statistic must be a function or the name of a built-in statistic.")
      }
    })
  } else {
    stop("statistics must be a character vector or a named list.")
  }

  # Step 2: transition counts, and for the parametric bootstrap the lengths
  # and first states of the sequences.
  if (is.matrix(data)) {
    if (method == "parametric") {
      stop("method = \"parametric\" needs the sequences, not only their counts.")
    }
    if (is.null(rownames(data)) || !identical(rownames(data), colnames(data)) ||
        any(!is.finite(data)) || any(data < 0) || any(data != round(data))) {
      stop("a count matrix must be square, with the states as identical row and ",
           "column names, and non-negative whole numbers.")
    }
    stateNames <- union(rownames(data), as.character(possibleStates))
    counts <- matrix(0, length(stateNames), length(stateNames),
                     dimnames = list(stateNames, stateNames))
    counts[rownames(data), colnames(data)] <- data
    sequences <- NULL
  } else {
    sequences <- if (is.list(data)) data else list(data)
    ok <- vapply(sequences, function(s) is.atomic(s) && !is.matrix(s) && length(s) >= 2L,
                 logical(1))
    if (!length(sequences) || !all(ok)) {
      stop("data must be a sequence, a list of sequences of length at least 2, ",
           "or a matrix of transition counts.")
    }
    sequences <- lapply(sequences, as.character)
    if (any(vapply(sequences, anyNA, logical(1)))) {
      stop("the sequences must not contain missing values.")
    }
    stateNames <- sort(unique(c(unlist(sequences), as.character(possibleStates))))
    counts <- Reduce(`+`, lapply(sequences, function(s) {
      createSequenceMatrix(s, possibleStates = stateNames)[stateNames, stateNames, drop = FALSE]
    }))
  }
  departures <- rowSums(counts)
  if (any(departures == 0)) {
    stop("the state(s) ", paste(stateNames[departures == 0], collapse = ", "),
         " are never left in the data, so their transitions cannot be estimated.")
  }

  # Step 3: the estimated chain and its statistics; their lengths fix the
  # columns of the replicate matrix.
  makeChain <- function(N) {
    new("markovchain", states = stateNames, byrow = TRUE,
        transitionMatrix = N / rowSums(N), name = "Bootstrap estimate")
  }
  chain <- makeChain(counts)
  estimates <- lapply(statistics, function(s) {
    value <- tryCatch(as.numeric(s(chain)), error = function(e) e)
    if (inherits(value, "error")) {
      stop("a statistic cannot be computed on the estimated chain: ", conditionMessage(value))
    }
    value
  })
  widths <- vapply(estimates, length, integer(1))
  if (any(widths == 0L)) {
    stop("every statistic must return at least one number.")
  }
  labels <- unlist(lapply(names(statistics), function(nm) {
    if (widths[[nm]] == 1L) {
      nm
    } else {
      own <- names(statistics[[nm]](chain))
      paste(nm, if (is.null(own)) seq_len(widths[[nm]]) else own, sep = ".")
    }
  }))
  estimate <- stats::setNames(unlist(estimates, use.names = FALSE), labels)

  # Step 4: the replicates.
  P <- chain@transitionMatrix
  replicateOnce <- function() {
    if (method == "multinomial") {
      N <- t(vapply(seq_along(stateNames), function(i)
        as.numeric(stats::rmultinom(1L, departures[i], P[i, ])), numeric(length(stateNames))))
      dimnames(N) <- dimnames(counts)
    } else {
      N <- Reduce(`+`, lapply(sequences, function(s) {
        path <- rmarkovchain(length(s) - 1L, chain, t0 = s[1L], include.t0 = TRUE)
        createSequenceMatrix(path, possibleStates = stateNames)[stateNames, stateNames, drop = FALSE]
      }))
      if (any(rowSums(N) == 0)) {
        return(rep(NA_real_, length(estimate)))
      }
    }
    .evaluateStatistics(makeChain(N), statistics, widths)
  }
  replicates <- .withLocalSeed(seed, t(vapply(seq_len(nboot), function(b) replicateOnce(),
                                              numeric(length(estimate)))))
  if (length(estimate) == 1L) {
    replicates <- matrix(replicates, ncol = 1L)
  }
  colnames(replicates) <- labels

  # Step 5: summaries and intervals.
  alpha <- 1 - level
  q <- apply(replicates, 2L, stats::quantile, probs = c(alpha / 2, 1 - alpha / 2),
             na.rm = TRUE, names = FALSE, type = 7)
  q <- matrix(q, nrow = 2L)
  ci <- if (type == "percentile") t(q) else cbind(2 * estimate - q[2L, ], 2 * estimate - q[1L, ])
  dimnames(ci) <- list(labels, c(sprintf("%g%%", 100 * alpha / 2),
                                 sprintf("%g%%", 100 * (1 - alpha / 2))))
  structure(list(
    estimate = estimate,
    ci = ci,
    se = apply(replicates, 2L, stats::sd, na.rm = TRUE),
    bias = colMeans(replicates, na.rm = TRUE) - estimate,
    failures = colSums(is.na(replicates)),
    replicates = replicates,
    chain = chain,
    counts = counts,
    method = method, type = type, level = level, nboot = nboot
  ), class = "markovchainBootstrap")
}

#' Print bootstrap intervals for a Markov chain
#'
#' @param x A \code{"markovchainBootstrap"} object returned by
#'   \code{\link{bootstrapProperties}}.
#' @param digits Number of significant digits.
#' @param ... Passed to \code{print}.
#' @return \code{x}, invisibly.
#' @method print markovchainBootstrap
#' @keywords internal
#' @export
print.markovchainBootstrap <- function(x, digits = 4, ...) {
  cat(sprintf("Bootstrap of a Markov chain on %d states (%s, %d replicates)\n\n",
              length(x$chain@states), x$method, x$nboot))
  table <- cbind(estimate = x$estimate, se = x$se, x$ci)
  print(signif(table, digits), ...)
  cat(sprintf("\n%s %g%% intervals.\n",
              if (x$type == "percentile") "Percentile" else "Basic", 100 * x$level))
  if (any(x$failures > 0)) {
    cat("Replicates where a statistic could not be computed (left out):",
        paste(sprintf("%s %d", names(x$failures)[x$failures > 0], x$failures[x$failures > 0]),
              collapse = ", "), "\n")
  }
  invisible(x)
}
