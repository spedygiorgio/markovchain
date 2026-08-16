#' Assess time-homogeneity of an empirical Markov chain
#'
#' @description
#' Tests whether the transition probabilities of a first-order Markov chain
#' are constant across a number of consecutive blocks. Structural zeros can
#' be supplied explicitly through a logical transition matrix. Sampling zeros
#' are retained as ordinary zero counts.
#'
#' @param sequence An empirical sequence.
#' @param nblocks Number of blocks.
#' @param structural.zeros Optional logical matrix indicating transitions that
#'   are structurally impossible. `TRUE` marks an impossible transition.
#'   Rows and columns correspond to the states in `sequence`. If dimnames are
#'   supplied, they are used to match states; otherwise the matrix is assumed
#'   to follow the order of `unique(sequence)`.
#' @param verbose Should test results be printed out?
#'
#' @return A list containing the chi-squared test statistic, degrees of freedom,
#'   and p-value.
#'
#' @references
#' Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about
#' Markov chains. The Annals of Mathematical Statistics, 28(1), 89-110.
#' @family statisticalTests
#' @rdname statisticalTests
#' @export
assessStationarity <- function(sequence, nblocks, structural.zeros = NULL,
                                verbose = TRUE) {
  warning("The accuracy of the statistical inference functions has been questioned. It will be thoroughly investigated in future versions of the package.")

  if (length(sequence) < 2L)
    stop("sequence must contain at least two observations")
  if (length(nblocks) != 1L || !is.numeric(nblocks) ||
      !is.finite(nblocks) || nblocks < 2 || nblocks != as.integer(nblocks))
    stop("nblocks must be a positive integer greater than one")

  nblocks <- as.integer(nblocks)
  states <- unique(sequence)
  nstates <- length(states)

  if (!is.null(structural.zeros)) {
    if (!is.matrix(structural.zeros) || !is.logical(structural.zeros))
      stop("structural.zeros must be a logical matrix")
    if (!all(dim(structural.zeros) == c(nstates, nstates)))
      stop("structural.zeros must have dimensions length(unique(sequence)) x length(unique(sequence))")
    if (anyNA(structural.zeros))
      stop("structural.zeros cannot contain NA values")

    if (!is.null(rownames(structural.zeros)) || !is.null(colnames(structural.zeros))) {
      if (is.null(rownames(structural.zeros)) || is.null(colnames(structural.zeros)) ||
          !setequal(rownames(structural.zeros), as.character(states)) ||
          !setequal(colnames(structural.zeros), as.character(states)))
        stop("rownames and colnames of structural.zeros must match the states in sequence")
      structural.zeros <- structural.zeros[as.character(states), as.character(states), drop = FALSE]
    }
  } else {
    structural.zeros <- matrix(FALSE, nstates, nstates)
  }

  # A transition marked structural zero must never occur in the data.
  observed <- matrix(FALSE, nstates, nstates,
                      dimnames = list(as.character(states), as.character(states)))
  for (j in seq_len(length(sequence) - 1L)) {
    from <- match(sequence[j], states)
    to <- match(sequence[j + 1L], states)
    observed[from, to] <- TRUE
  }
  if (any(observed & structural.zeros))
    stop("sequence contains a transition marked as a structural zero")

  blocksize <- length(sequence) / nblocks
  TStat <- 0
  df <- 0

  for (i in seq_len(nstates)) {
    counts <- matrix(0, nrow = nblocks, ncol = nstates,
                     dimnames = list(seq_len(nblocks), as.character(states)))

    for (j in seq_len(length(sequence) - 1L)) {
      if (sequence[j] == states[i]) {
        block <- min(nblocks, ceiling(j / blocksize))
        future <- match(sequence[j + 1L], states)
        counts[block, future] <- counts[block, future] + 1
      }
    }

    # A block with no observations from state i carries no information about
    # the conditional distribution P(X[t+1] | X[t] = i).
    counts <- counts[rowSums(counts) > 0, , drop = FALSE]

    # Remove only transitions explicitly declared structurally impossible.
    # Sampling zeros remain in the table and therefore contribute correctly
    # to the Pearson statistic and its degrees of freedom.
    counts <- counts[, !structural.zeros[i, ], drop = FALSE]

    if (nrow(counts) < 2L || ncol(counts) < 2L)
      next

    expected <- outer(rowSums(counts), colSums(counts) / sum(counts))
    positive.expected <- expected > 0
    statistic.i <- sum((counts[positive.expected] - expected[positive.expected])^2 /
                       expected[positive.expected])
    df.i <- (nrow(counts) - 1L) * (ncol(counts) - 1L)

    TStat <- TStat + statistic.i
    df <- df + df.i
  }

  if (df <= 0) {
    pvalue <- NA_real_
    warning("Insufficient non-structural observations to perform the stationarity test")
  } else {
    pvalue <- pchisq(TStat, df = df, lower.tail = FALSE)
  }

  if (verbose) {
    cat("The assessStationarity test statistic is: ", TStat, "\n")
    cat("The Chi-Square d.f. are: ", df, "\n")
    cat("The p-value is: ", pvalue, "\n")
  }

  invisible(list(statistic = TStat, dof = df, p.value = pvalue))
}
