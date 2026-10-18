# Plots of state sequences: a sequence index plot (one row per sequence, one
# coloured cell per time point, the "heatmap of trajectories") and a state
# distribution plot (the proportion of sequences in each state at each time,
# also called a chronogram), with the distribution expected from a chain
# drawn on top.
#
# For a junior reader: a transition matrix summarises what happens in one
# step; these plots show whole trajectories, so that one can see at a glance
# how long the chain stays in a state, whether some trajectories behave
# differently, and whether the observed state distribution evolves over time
# as the fitted chain predicts.

# Internal helper: turn the data argument into a character matrix with one
# row per sequence and one column per time point (shorter sequences padded
# with NA), and the ordered vector of states.
.sequenceMatrix <- function(x, stateOrder = NULL) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (is.matrix(x)) {
    levelsFromFactor <- NULL
    m <- matrix(as.character(x), nrow(x), ncol(x))
  } else {
    sequences <- if (is.list(x)) x else list(x)
    ok <- vapply(sequences, function(s) is.atomic(s) && length(s) >= 1L, logical(1))
    if (!length(sequences) || !all(ok)) {
      stop("x must be a markovchain, a sequence, a list of sequences, or a matrix or ",
           "data frame with one sequence per row.")
    }
    levelsFromFactor <- unique(unlist(lapply(sequences, function(s) if (is.factor(s)) levels(s))))
    len <- max(lengths(sequences))
    m <- t(vapply(sequences, function(s) {
      s <- as.character(s)
      c(s, rep(NA_character_, len - length(s)))
    }, character(len)))
    if (len == 1L) {
      m <- matrix(m, ncol = 1L)
    }
  }
  observed <- unique(m[!is.na(m)])
  if (!length(observed)) {
    stop("the sequences contain no states.")
  }
  if (!is.null(stateOrder)) {
    unknown <- setdiff(observed, stateOrder)
    if (length(unknown)) {
      stop("states not in the chain: ", paste(unknown, collapse = ", "), ".")
    }
    stateNames <- stateOrder
  } else if (length(levelsFromFactor)) {
    stateNames <- union(levelsFromFactor, sort(observed))
  } else {
    stateNames <- sort(observed)
  }
  list(sequences = m, states = stateNames)
}

#' Plot state sequences of a Markov chain
#'
#' Draws simulated or observed trajectories of a Markov chain either as a
#' sequence index plot (each sequence a row of coloured cells, one per time
#' point) or as a state distribution plot (the proportion of sequences in
#' each state at each time point), optionally with the distribution expected
#' from a chain overlaid.
#'
#' @param x A \code{markovchain} object, from which \code{nSequences}
#'   sequences of \code{steps} time points are simulated, or observed data: a
#'   sequence, a list of sequences (of possibly different lengths), or a
#'   matrix or data frame with one sequence per row and one time point per
#'   column. Missing values are left blank.
#' @param type \code{"index"} (the default) for a sequence index plot,
#'   \code{"distribution"} for a state distribution plot.
#' @param nSequences,steps For a \code{markovchain}: the number of sequences
#'   to simulate and their number of time points (the initial state
#'   included).
#' @param initial For a \code{markovchain}: the distribution of the initial
#'   state, as in \code{\link{redistribute}}: \code{NULL} (uniform, the
#'   default), a state name, or a numeric vector of probabilities.
#' @param chain For observed data, an optional \code{markovchain} on the same
#'   states whose expected state distribution, started from the observed
#'   distribution at the first time point, is overlaid on the state
#'   distribution plot. For a \code{markovchain} \code{x} it is \code{x}
#'   itself.
#' @param sort For \code{type = "index"}: \code{"none"} (the default) keeps
#'   the sequences in their order; \code{"first"} and \code{"last"} sort them
#'   by their state at the first or last time point (ties broken by the
#'   following or preceding time points).
#' @param expected For \code{type = "distribution"}: whether to overlay the
#'   expected distribution when a chain is available (default \code{TRUE}).
#' @param col Optional vector of colours, one per state. The default is
#'   \code{grDevices::hcl.colors(n, "Dark 3")}.
#' @param legend Whether to draw a legend of the states (default
#'   \code{TRUE}).
#' @param main,xlab,ylab Title and axis labels.
#' @param seed For a \code{markovchain}: an optional whole number. When given,
#'   the sequences are simulated after \code{set.seed(seed)} and the caller's
#'   random number stream is left untouched.
#' @param ... Further graphical parameters passed to \code{plot.default}.
#'
#' @return Invisibly, a list with \code{sequences}, the character matrix of
#'   the sequences as plotted (one row per sequence, after sorting);
#'   \code{proportions}, the matrix of the observed proportion of each state
#'   (rows) at each time point (columns); and \code{expected}, the matching
#'   matrix of expected proportions, or \code{NULL} when no chain is
#'   available.
#'
#' @details
#' Sequence index plots and state distribution plots are the two basic
#' graphs of sequence analysis (Gabadinho et al. 2011). For a Markov chain
#' the first shows the lengths of the spells in each state and the second
#' how the distribution of the state moves towards the stationary
#' distribution. The expected distribution at time \eqn{t} is
#' \eqn{\mu_1 P^{t-1}}, with \eqn{\mu_1} the observed (or, for a
#' \code{markovchain} \code{x}, the given) initial distribution, computed by
#' \code{\link{redistribute}}; with observed data a systematic gap between the
#' two suggests that the chain is not time-homogeneous (see
#' \code{\link{assessStationarity}}). The proportions at each time point
#' are computed over the sequences observed at that time point.
#'
#' The plot is drawn with base graphics. The legend is placed in the right
#' margin, which is widened for the purpose and restored afterwards.
#'
#' @references
#' Gabadinho, A., Ritschard, G., Mueller, N. S. and Studer, M. (2011).
#' Analyzing and visualizing state sequences in R with TraMineR.
#' \emph{Journal of Statistical Software}, 40(4), 1-37.
#' \doi{10.18637/jss.v040.i04}
#'
#' @seealso \code{\link{rmarkovchain}}, \code{\link{redistribute}},
#'   \code{\link{autoplot.markovchain}} (\code{type = "flow"})
#'
#' @examples
#' s <- c("sunny", "cloudy", "rain")
#' weather <- new("markovchain", states = s, transitionMatrix = matrix(
#'   c(0.7, 0.2, 0.1, 0.3, 0.4, 0.3, 0.2, 0.45, 0.35), 3, byrow = TRUE,
#'   dimnames = list(s, s)))
#' plotSequence(weather, nSequences = 30, steps = 40, initial = "rain", seed = 1)
#' plotSequence(weather, type = "distribution", nSequences = 200, steps = 15,
#'              initial = "rain", seed = 1)
#'
#' # observed trajectories, sorted by their first state, and their
#' # distribution against a time-homogeneous chain fitted to them
#' data(holson)
#' plotSequence(holson[1:60, 2:12], sort = "first")
#' fitted <- markovchainFit(holson[, 2:12])$estimate
#' plotSequence(holson[, 2:12], type = "distribution", chain = fitted)
#'
#' @export
plotSequence <- function(x, type = c("index", "distribution"), nSequences = 20L, steps = 50L,
                         initial = NULL, chain = NULL, sort = c("none", "first", "last"),
                         expected = TRUE, col = NULL, legend = TRUE, main = NULL,
                         xlab = "Time", ylab = NULL, seed = NULL, ...) {
  type <- match.arg(type)
  sort <- match.arg(sort)

  # Step 1: the sequences, simulated or observed.
  if (is(x, "markovchain")) {
    for (arg in c("nSequences", "steps")) {
      v <- get(arg)
      if (length(v) != 1L || !is.numeric(v) || is.na(v) || v != as.integer(v) || v < 1L) {
        stop(arg, " must be a single positive integer.")
      }
    }
    stateNames <- x@states
    mu <- .initialDistribution(initial, stateNames)
    P <- .rowStochasticMatrix(x)
    m <- .withLocalSeed(seed, t(vapply(seq_len(nSequences), function(k) {
      out <- character(steps)
      cur <- sample.int(length(stateNames), 1L, prob = mu)
      out[1L] <- stateNames[cur]
      for (t in seq_len(steps - 1L)) {
        cur <- sample.int(length(stateNames), 1L, prob = P[cur, ])
        out[t + 1L] <- stateNames[cur]
      }
      out
    }, character(steps))))
    if (steps == 1L) {
      m <- matrix(m, ncol = 1L)
    }
    chain <- x
    if (is.null(main)) {
      main <- if (type == "index") paste("Simulated sequences of", x@name)
              else paste("State distribution of", x@name)
    }
  } else {
    if (!is.null(chain) && !is(chain, "markovchain")) {
      stop("chain must be a markovchain object.")
    }
    sm <- .sequenceMatrix(x, if (is.null(chain)) NULL else chain@states)
    m <- sm$sequences
    stateNames <- sm$states
    if (is.null(main)) {
      main <- if (type == "index") "Sequence index plot" else "State distribution plot"
    }
  }
  k <- length(stateNames)
  nTimes <- ncol(m)
  if (is.null(col)) {
    col <- grDevices::hcl.colors(max(k, 2L), "Dark 3")[seq_len(k)]
  } else if (length(col) < k) {
    stop("col must give at least one colour per state.")
  }
  codes <- matrix(match(m, stateNames), nrow(m), nTimes)

  # Step 2: observed proportions and, if a chain is available, expected ones.
  counts <- vapply(seq_len(nTimes), function(t) tabulate(codes[, t], nbins = k), numeric(k))
  counts <- matrix(counts, k, nTimes, dimnames = list(stateNames, NULL))
  totals <- colSums(counts)
  proportions <- sweep(counts, 2L, pmax(totals, 1), "/")
  proportions[, totals == 0] <- NA
  expectedProportions <- NULL
  if (!is.null(chain)) {
    if (!setequal(chain@states, stateNames)) {
      stop("chain must have the same states as the sequences.")
    }
    start <- if (is(x, "markovchain")) {
      stats::setNames(mu, stateNames)
    } else {
      if (totals[1L] == 0) stop("no sequence is observed at the first time point.")
      stats::setNames(proportions[, 1L], stateNames)
    }
    traj <- redistribute(chain, steps = nTimes - 1L, initial = start[chain@states])
    traj <- matrix(traj, ncol = length(chain@states), dimnames = list(NULL, chain@states))
    expectedProportions <- t(traj[, stateNames, drop = FALSE])
  }

  # Step 3: sorting of the sequences for the index plot.
  if (sort != "none") {
    keys <- as.data.frame(codes)
    if (sort == "last") {
      keys <- keys[, rev(seq_len(nTimes)), drop = FALSE]
    }
    ord <- do.call(order, c(unname(as.list(keys)), list(na.last = TRUE)))
    m <- m[ord, , drop = FALSE]
    codes <- codes[ord, , drop = FALSE]
  }

  # Step 4: draw. The right margin holds the legend.
  oldPar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(mar = oldPar$mar, xpd = oldPar$xpd))
  if (legend) {
    # widen the right margin to the longest label plus the colour box, in
    # lines of text
    labels <- c(stateNames, if (type == "distribution" && expected &&
                                !is.null(expectedProportions)) "expected")
    needed <- (max(graphics::strwidth(labels, units = "inches", cex = 0.85)) + 0.45) /
      graphics::par("csi") + 0.5
    mar <- oldPar$mar
    mar[4L] <- max(mar[4L], needed)
    graphics::par(mar = mar)
  }
  times <- seq_len(nTimes)
  if (type == "index") {
    nSeq <- nrow(codes)
    if (is.null(ylab)) ylab <- "Sequence"
    graphics::plot.default(NA, xlim = c(0.5, nTimes + 0.5), ylim = c(0.5, nSeq + 0.5),
                           xlab = xlab, ylab = ylab, main = main, yaxs = "i", xaxs = "i",
                           ...)
    # the first sequence at the top, as in a table
    graphics::image(times, seq_len(nSeq), t(codes[rev(seq_len(nSeq)), , drop = FALSE]),
                    col = col[seq_len(k)], breaks = seq(0.5, k + 0.5), add = TRUE,
                    useRaster = nSeq * nTimes > 5000)
    graphics::box()
  } else {
    if (is.null(ylab)) ylab <- "Proportion"
    graphics::plot.default(NA, xlim = range(times) + c(-0.5, 0.5), ylim = c(0, 1),
                           xlab = xlab, ylab = ylab, main = main, xaxs = "i", yaxs = "i", ...)
    cumObs <- apply(rbind(0, proportions), 2L, cumsum)
    cumObs <- matrix(cumObs, k + 1L, nTimes)
    for (s in seq_len(k)) {
      for (t in times[!is.na(cumObs[1L + s, ])]) {
        graphics::rect(t - 0.5, cumObs[s, t], t + 0.5, cumObs[s + 1L, t],
                       col = col[s], border = NA)
      }
    }
    if (expected && !is.null(expectedProportions)) {
      cumExp <- apply(expectedProportions, 2L, cumsum)
      cumExp <- matrix(cumExp, k, nTimes)
      for (s in seq_len(k - 1L)) {
        graphics::lines(times, cumExp[s, ], lty = 2, lwd = 2)
      }
    }
    graphics::box()
  }
  if (legend) {
    graphics::par(xpd = TRUE)
    usr <- graphics::par("usr")
    graphics::legend(usr[2L] + 0.02 * diff(usr[1:2]), usr[4L], legend = stateNames,
                     fill = col[seq_len(k)], bty = "n", cex = 0.85)
    if (type == "distribution" && expected && !is.null(expectedProportions)) {
      graphics::legend(usr[2L] + 0.02 * diff(usr[1:2]), usr[3L] + 0.2 * diff(usr[3:4]),
                       legend = "expected", lty = 2, lwd = 2, bty = "n", cex = 0.85)
    }
  }
  invisible(list(sequences = m, proportions = proportions, expected = expectedProportions))
}
