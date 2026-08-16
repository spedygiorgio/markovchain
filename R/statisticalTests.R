# Revised statistical inference functions for empirical Markov chains.

# This file follows the same high-level pattern for all tests:
# validate inputs -> build observed counts -> compute a test statistic and
# degrees of freedom -> optionally simulate the null distribution -> return
# a standard `htest` object.  The helper below centralizes the final formatting.

#' Standardize statistical inference results to R's htest convention.
#'
#' @keywords internal
.asHtest <- function(result, method, data.name) {
  # Convert the package-specific result into the structure expected by R's
  # standard hypothesis-testing methods.
  statistic.name <- switch(method, Pearson = "X-squared", G = "G-squared",
                            simulation = "G-squared")
  method.name <- switch(method, Pearson = "Pearson's Chi-squared test",
                        G = "Likelihood-ratio test",
                        simulation = "Parametric Monte Carlo test (likelihood-ratio statistic)")
  result$statistic <- stats::setNames(as.numeric(result$statistic), statistic.name)
  result$parameter <- stats::setNames(as.numeric(result$dof), "df")
  result$method <- method.name
  result$data.name <- data.name
  class(result) <- c("htest", setdiff(class(result), "htest"))
  result
}

#' Test the first-order Markov property of an empirical sequence
#'
#' Tests the null hypothesis that the conditional distribution of the next
#' state depends only on the current state, against a second-order alternative.
#' @rdname statisticalTests
#' @family statisticalTests
#' @param sequence An empirical sequence of states.
#' @param method Test statistic: `"G"`, `"Pearson"`, or `"simulation"`.
#' @param B Number of Monte Carlo replicates for simulation.
#' @param seed Optional random seed.
#' @param verbose Should results be printed?
#' @return An `htest` object with additional package-specific components.
#' @references Anderson, T. W. and Goodman, L. A. (1957). Statistical inference
#' about Markov chains. *The Annals of Mathematical Statistics*, 28(1), 89--110.
#' @export
verifyMarkovProperty <- function(sequence, method = c("G", "Pearson", "simulation"),
                                 B = 9999, seed = NULL, verbose = TRUE) {
  method <- match.arg(method)
  if (length(sequence) < 4L) stop("sequence must contain at least four observations.")
  if (anyNA(sequence)) stop("sequence must not contain missing values.")
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B))
    stop("B must be a positive integer.")

  # Encode the states as integer indices.  This makes the three-way transition
  # table N(i,j,k) compact and allows it to be populated with `tabulate()`.
  states <- unique(sequence); r <- length(states); n <- length(sequence)
  idx <- match(sequence, states)
  Nijk <- array(0L, dim = c(r, r, r), dimnames = list(states, states, states))
  linear <- idx[-c(1L, 2L)] + (idx[2:(n - 1L)] - 1L) * r +
    (idx[1:(n - 2L)] - 1L) * r * r
  Nijk[] <- tabulate(linear, nbins = r^3L)

  # Collapse the three-way table to the two margins needed by the null model
  # and construct its expected counts.
  Nij <- apply(Nijk, c(1L, 2L), sum); Njk <- apply(Nijk, c(2L, 3L), sum)
  Nj <- rowSums(Njk)
  expected <- array(0, dim = c(r, r, r), dimnames = dimnames(Nijk))
  for (j in seq_len(r)) if (Nj[j] > 0L)
    expected[, j, ] <- outer(Nij[, j], Njk[j, ]) / Nj[j]

  # Use one statistic implementation for both the Pearson and likelihood-ratio
  # versions, taking care not to divide by or take logs of zero expected counts.
  statistic_fun <- function(observed, expected, method) {
    positive_expected <- expected > 0
    if (method == "Pearson")
      return(sum((observed[positive_expected] - expected[positive_expected])^2 /
                 expected[positive_expected]))
    positive <- positive_expected & observed > 0
    2 * sum(observed[positive] * log(observed[positive] / expected[positive]))
  }
  statistic <- statistic_fun(Nijk, expected, if (method == "Pearson") "Pearson" else "G")

  # The degrees of freedom are obtained by summing the contingency-table
  # degrees of freedom over each intermediate state j.
  dof <- 0L
  for (j in seq_len(r)) {
    previous <- sum(Nij[, j] > 0); future <- sum(Njk[j, ] > 0)
    if (previous > 0L && future > 0L) dof <- dof + (previous - 1L) * (future - 1L)
  }

  # Estimate the first-order transition matrix.  It is also the transition
  # mechanism used to generate sequences under the null in simulation mode.
  P <- matrix(0, r, r, dimnames = list(states, states)); informative <- Nj > 0
  P[informative, ] <- sweep(Njk[informative, , drop = FALSE], 1L, Nj[informative], "/")
  if (method != "simulation") {
    # Asymptotic chi-squared p-value.
    p.value <- if (dof > 0) pchisq(statistic, dof, lower.tail = FALSE) else NA_real_
    result <- list(statistic = statistic, dof = dof, p.value = p.value, method = method,
                   observed = Nijk, expected = expected, transitionMatrix = P)
  } else {
    # Simulation requires every state in the sequence to have an observed
    # outgoing transition, so that P defines a valid sampling distribution.
    if (!all(informative)) stop("method = 'simulation' requires at least one observed outgoing transition from every state appearing in the sequence.")
    if (!is.null(seed)) set.seed(seed)
    simulated <- numeric(B)
    for (b in seq_len(B)) {
      # Generate a sequence from the fitted first-order Markov chain.
      x <- integer(n); x[1L] <- idx[1L]
      for (t in seq_len(n - 1L)) x[t + 1L] <- sample.int(r, 1L, prob = P[x[t], ])
      a <- x[-c(1L, 2L)] + (x[2:(n - 1L)] - 1L) * r + (x[1:(n - 2L)] - 1L) * r * r
      Ns <- array(tabulate(a, nbins = r^3L), dim = c(r, r, r))
      Nijs <- apply(Ns, c(1L, 2L), sum); Njks <- apply(Ns, c(2L, 3L), sum); Njs <- rowSums(Njks)
      Es <- array(0, dim = c(r, r, r))
      for (j in seq_len(r)) if (Njs[j] > 0L) Es[, j, ] <- outer(Nijs[, j], Njks[j, ]) / Njs[j]
      simulated[b] <- statistic_fun(Ns, Es, if (method == "Pearson") "Pearson" else "G")
    }
    # Add-one correction gives the Monte Carlo p-value and avoids returning
    # exactly zero when no simulated statistic exceeds the observed one.
    p.value <- (1 + sum(simulated >= statistic)) / (B + 1)
    result <- list(statistic = statistic, dof = dof, p.value = p.value, method = method,
                   observed = Nijk, expected = expected, transitionMatrix = P,
                   simulations = simulated, B = B, initial.state = states[idx[1L]])
  }
  result <- .asHtest(result, method, deparse(substitute(sequence)))
  if (verbose) print(result)
  invisible(result)
}

#' Assess the order of an empirical Markov chain
#' @rdname statisticalTests
#' @family statisticalTests
#' @param sequence An empirical sequence of states.
#' @param verbose Should test results be printed?
#' @return An `htest` object.
#' @export
assessOrder <- function(sequence, verbose = TRUE) {
  method <- "Pearson"
  if (length(sequence) < 4L) stop("sequence must contain at least four observations.")
  if (anyNA(sequence)) stop("sequence must not contain missing values.")
  states <- unique(sequence); k <- length(states); n <- length(sequence); statistic <- 0

  # For each possible intermediate state, compare the observed two-step
  # transitions with the distribution expected under a first-order chain.
  for (present in states) {
    mat <- matrix(0, nrow = k, ncol = k, dimnames = list(states, states))
    for (i in seq_len(n - 2L)) if (identical(sequence[i + 1L], present)) {
      mat[as.character(sequence[i]), as.character(sequence[i + 2L])] <-
        mat[as.character(sequence[i]), as.character(sequence[i + 2L])] + 1
    }
    row_totals <- rowSums(mat); active <- row_totals > 0
    if (sum(active) > 1L && sum(colSums(mat) > 0) > 1L) {
      expected <- outer(row_totals[active], colSums(mat)) / sum(row_totals[active])
      observed <- mat[active, , drop = FALSE]; positive <- expected > 0
      statistic <- statistic + sum((observed[positive] - expected[positive])^2 / expected[positive])
    }
  }
  dof <- k * (k - 1L)^2
  result <- list(statistic = statistic, dof = dof,
                 p.value = if (dof > 0L) pchisq(statistic, dof, lower.tail = FALSE) else NA_real_)
  result <- .asHtest(result, method, deparse(substitute(sequence)))
  result$method <- "Pearson's Chi-squared test for Markov order"
  if (verbose) print(result)
  invisible(result)
}

#' Test compatibility of an empirical transition matrix with a theoretical one
#' @rdname statisticalTests
#' @family statisticalTests
#' @param data An empirical sequence or a matrix of transition counts.
#' @param object A `markovchain` object specifying theoretical probabilities.
#' @param method Test statistic: `"G"`, `"Pearson"`, or `"simulation"`.
#' @param B Number of Monte Carlo replicates for simulation.
#' @param seed Optional random seed.
#' @param verbose Should results be printed?
#' @return An `htest` object with observed and expected counts.
#' @references Kullback, S., Kupperman, M. and Ku, H. H. (1962). Tests for Contingency Tables and Markov Chains. *Technometrics*, 4(4), 573--608.
#' @export
verifyEmpiricalToTheoretical <- function(data, object, method = c("G", "Pearson", "simulation"),
                                         B = 9999, seed = NULL, verbose = TRUE) {
  method <- match.arg(method)
  if (missing(data) || missing(object)) stop("Both data and object must be supplied.")
  if (!methods::is(object, "markovchain")) stop("object must be a markovchain object.")
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B)) stop("B must be a positive integer.")
  states <- names(object); P <- object@transitionMatrix

  # Accept either a transition-count matrix or a raw empirical sequence and
  # convert both representations to a common count matrix.
  if (is.matrix(data)) {
    observed <- data
    if (is.null(rownames(observed)) || is.null(colnames(observed))) stop("A transition-count matrix must have row and column names.")
  } else if (is.atomic(data)) observed <- createSequenceMatrix(stringchar = data, possibleStates = states)
  else stop("data must be an empirical sequence or transition-count matrix.")
  if (anyNA(observed) || any(observed < 0) || any(observed != floor(observed))) stop("Observed transition counts must be non-negative integers.")
  if (!setequal(rownames(observed), states) || !setequal(colnames(observed), states)) stop("Empirical and theoretical transition matrices must have the same state support.")

  # Reorder the observed matrix to match the theoretical state order, then
  # calculate expected counts from the theoretical transition probabilities.
  observed <- observed[states, states, drop = FALSE]; row_totals <- rowSums(observed)
  expected <- sweep(P, 1L, row_totals, "*"); informative <- row_totals > 0
  if (any(observed > 0 & P == 0)) {
    # An observed transition that is theoretically impossible gives an
    # infinite likelihood-ratio statistic and a p-value of zero.
    result <- list(statistic = Inf, dof = NA_integer_, p.value = 0, method = method,
                   observed = observed, expected = expected, transitionMatrix = P)
  } else {
    statistic <- if (method == "Pearson") {
      pos <- outer(informative, rep(TRUE, length(states))) & P > 0
      sum((observed[pos] - expected[pos])^2 / expected[pos])
    } else {
      pos <- P > 0 & observed > 0
      2 * sum(observed[pos] * log(observed[pos] / expected[pos]))
    }
    dof <- sum(vapply(which(informative), function(i) sum(P[i, ] > 0) - 1L, integer(1)))
    if (method != "simulation") {
      p.value <- if (dof > 0) pchisq(statistic, dof, lower.tail = FALSE) else NA_real_
      result <- list(statistic = statistic, dof = dof, p.value = p.value, method = method,
                     observed = observed, expected = expected, transitionMatrix = P)
    } else {
      if (!is.null(seed)) set.seed(seed); simulated <- numeric(B)
      # Conditional on each row total, simulate multinomial counts using the
      # theoretical transition probabilities as the null distribution.
      for (b in seq_len(B)) {
        sim <- matrix(0, nrow(P), ncol(P), dimnames = dimnames(P))
        for (i in which(informative)) sim[i, ] <- as.vector(rmultinom(1L, row_totals[i], P[i, ]))
        if (method == "Pearson") {
          pos <- outer(informative, rep(TRUE, length(states))) & P > 0
          simulated[b] <- sum((sim[pos] - expected[pos])^2 / expected[pos])
        } else {
          pos <- P > 0 & sim > 0
          simulated[b] <- 2 * sum(sim[pos] * log(sim[pos] / expected[pos]))
        }
      }
      p.value <- (1 + sum(simulated >= statistic)) / (B + 1)
      result <- list(statistic = statistic, dof = dof, p.value = p.value, method = method,
                     observed = observed, expected = expected, transitionMatrix = P,
                     simulations = simulated, B = B)
    }
  }
  result <- .asHtest(result, method, deparse(substitute(data)))
  if (verbose) print(result)
  invisible(result)
}

#' Test homogeneity of transition probabilities across empirical sequences
#' @rdname statisticalTests
#' @family statisticalTests
#' @param inputList A list whose elements are empirical sequences or matrices.
#' @param method Test statistic: `"G"`, `"Pearson"`, or `"simulation"`.
#' @param B Number of Monte Carlo replicates for simulation.
#' @param seed Optional random seed.
#' @param verbose Should results be printed?
#' @return An `htest` object with pooled and individual transition counts.
#' @export
verifyHomogeneity <- function(inputList, method = c("G", "Pearson", "simulation"),
                              B = 9999, seed = NULL, verbose = TRUE) {
  method <- match.arg(method)
  if (!is.list(inputList) || length(inputList) < 2L) stop("inputList must contain at least two sequences or matrices.")
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B)) stop("B must be a positive integer.")

  # Normalize each input to a transition-count matrix.  Missing states are
  # added as zero rows/columns below so all samples share the same support.
  mats <- lapply(inputList, function(x) {
    if (is.matrix(x)) {
      if (is.null(rownames(x)) || is.null(colnames(x))) stop("Transition-count matrices must have row and column names.")
      m <- x
    } else if (is.atomic(x)) m <- createSequenceMatrix(stringchar = x)
    else stop("Each element must be a sequence or transition-count matrix.")
    if (anyNA(m) || any(m < 0) || any(m != floor(m))) stop("Transition counts must be non-negative integers.")
    m
  })
  states <- sort(unique(unlist(lapply(mats, rownames))))
  mats <- lapply(mats, function(m) {
    out <- matrix(0, length(states), length(states), dimnames = list(states, states)); out[rownames(m), colnames(m)] <- m; out
  })
  pooled <- Reduce("+", mats)

  # Compute the homogeneity statistic row by row.  For each origin state,
  # the null hypothesis is that all input samples share the same destination
  # probabilities; samples with no observations from that state are ignored.
  hom_stat <- function(mats0) {
    pooled0 <- Reduce("+", mats0); stat <- 0
    for (i in seq_along(states)) {
      totals <- vapply(mats0, function(m) sum(m[i, ]), numeric(1)); active <- totals > 0
      if (sum(active) < 2L) next
      destinations <- pooled0[i, ] > 0; k <- sum(destinations); if (k < 2L) next
      p <- pooled0[i, destinations] / sum(pooled0[i, destinations])
      for (s in which(active)) {
        obs <- mats0[[s]][i, destinations]; exp <- totals[s] * p
        if (method == "Pearson") stat <- stat + sum((obs - exp)^2 / exp)
        else { pos <- obs > 0; stat <- stat + 2 * sum(obs[pos] * log(obs[pos] / exp[pos])) }
      }
    }
    stat
  }
  statistic <- hom_stat(mats); dof <- 0L
  for (i in seq_along(states)) {
    totals <- vapply(mats, function(m) sum(m[i, ]), numeric(1)); active <- totals > 0
    destinations <- pooled[i, ] > 0; k <- sum(destinations)
    if (sum(active) >= 2L && k >= 2L) dof <- dof + (sum(active) - 1L) * (k - 1L)
  }
  if (method != "simulation") {
    p.value <- if (dof > 0) pchisq(statistic, dof, lower.tail = FALSE) else NA_real_
  } else {
    if (!is.null(seed)) set.seed(seed); simulated <- numeric(B)
    for (b in seq_len(B)) {
      # Simulate each sample independently, conditional on its observed row
      # totals, but using the pooled transition probabilities under H0.
      sim_mats <- lapply(mats, function(m) matrix(0, nrow(m), ncol(m), dimnames = dimnames(m)))
      for (i in seq_along(states)) {
        destinations <- pooled[i, ] > 0; if (!any(destinations)) next
        p <- pooled[i, destinations] / sum(pooled[i, destinations])
        for (s in seq_along(mats)) {
          total_s <- sum(mats[[s]][i, ])
          if (total_s > 0) sim_mats[[s]][i, destinations] <- as.vector(rmultinom(1L, total_s, p))
        }
      }
      simulated[b] <- hom_stat(sim_mats)
    }
    p.value <- (1 + sum(simulated >= statistic)) / (B + 1)
  }
  result <- list(statistic = statistic, dof = dof, p.value = p.value, method = method,
                 pooled = pooled, transitionCounts = mats)
  if (method == "simulation") { result$simulations <- simulated; result$B <- B }
  result <- .asHtest(result, method, deparse(substitute(inputList)))
  if (verbose) print(result)
  invisible(result)
}

#' Assess time-homogeneity of an empirical Markov chain
#'
#' Tests whether transition probabilities are constant across consecutive blocks.
#' Structural zeros can be supplied explicitly through a logical transition matrix.
#' @rdname statisticalTests
#' @family statisticalTests
#' @param sequence An empirical sequence.
#' @param nblocks Number of blocks, at least two.
#' @param structural.zeros Optional logical matrix marking impossible transitions.
#' @param verbose Should test results be printed?
#' @return An `htest` object.
#' @references Anderson, T. W. and Goodman, L. A. (1957). Statistical inference
#' about Markov chains. *The Annals of Mathematical Statistics*, 28(1), 89--110.
#' @export
assessStationarity <- function(sequence, nblocks, structural.zeros = NULL, verbose = TRUE) {
  if (length(sequence) < 2L) stop("sequence must contain at least two observations")
  if (length(nblocks) != 1L || !is.numeric(nblocks) || !is.finite(nblocks) ||
      nblocks < 2 || nblocks != as.integer(nblocks))
    stop("nblocks must be a positive integer greater than one")
  nblocks <- as.integer(nblocks); states <- unique(sequence); nstates <- length(states)

  # `structural.zeros` distinguishes impossible transitions from transitions
  # that simply have zero observed frequency.  Such cells are removed before
  # calculating the chi-squared statistic, preventing artificial zero expected
  # counts from producing NaN values.
  if (!is.null(structural.zeros)) {
    if (!is.matrix(structural.zeros) || !is.logical(structural.zeros)) stop("structural.zeros must be a logical matrix")
    if (!all(dim(structural.zeros) == c(nstates, nstates))) stop("structural.zeros must have dimensions length(unique(sequence)) x length(unique(sequence))")
    if (anyNA(structural.zeros)) stop("structural.zeros cannot contain NA values")
    if (!is.null(rownames(structural.zeros)) || !is.null(colnames(structural.zeros))) {
      if (is.null(rownames(structural.zeros)) || is.null(colnames(structural.zeros)) ||
          !setequal(rownames(structural.zeros), as.character(states)) || !setequal(colnames(structural.zeros), as.character(states)))
        stop("rownames and colnames of structural.zeros must match the states in sequence")
      structural.zeros <- structural.zeros[as.character(states), as.character(states), drop = FALSE]
    }
  } else structural.zeros <- matrix(FALSE, nstates, nstates)

  # Record which transitions are observed at least once.  This is used to
  # reject data that contradict the structural-zero specification.
  observed <- matrix(FALSE, nstates, nstates, dimnames = list(as.character(states), as.character(states)))
  for (j in seq_len(length(sequence) - 1L)) {
    from <- match(sequence[j], states); to <- match(sequence[j + 1L], states); observed[from, to] <- TRUE
  }
  if (any(observed & structural.zeros)) stop("sequence contains a transition marked as a structural zero")

  # Split the sequence into consecutive blocks.  For each origin state we
  # obtain a contingency table whose rows are time blocks and whose columns
  # are admissible destination states.
  blocksize <- length(sequence) / nblocks; TStat <- 0; df <- 0
  for (i in seq_len(nstates)) {
    counts <- matrix(0, nrow = nblocks, ncol = nstates, dimnames = list(seq_len(nblocks), as.character(states)))
    for (j in seq_len(length(sequence) - 1L)) if (sequence[j] == states[i]) {
      block <- min(nblocks, ceiling(j / blocksize)); future <- match(sequence[j + 1L], states); counts[block, future] <- counts[block, future] + 1
    }
    counts <- counts[rowSums(counts) > 0, , drop = FALSE]; counts <- counts[, !structural.zeros[i, ], drop = FALSE]
    if (nrow(counts) < 2L || ncol(counts) < 2L) next

    # Under time-homogeneity, the destination distribution is the same in
    # every non-empty block.  Expected counts are therefore block totals times
    # the pooled destination proportions.
    expected <- outer(rowSums(counts), colSums(counts) / sum(counts)); positive.expected <- expected > 0
    TStat <- TStat + sum((counts[positive.expected] - expected[positive.expected])^2 / expected[positive.expected])
    df <- df + (nrow(counts) - 1L) * (ncol(counts) - 1L)
  }
  pvalue <- if (df > 0) pchisq(TStat, df = df, lower.tail = FALSE) else NA_real_
  result <- list(statistic = TStat, dof = df, p.value = pvalue)
  result <- .asHtest(result, "Pearson", deparse(substitute(sequence)))
  result$method <- "Pearson's Chi-squared test for time-homogeneity"
  if (verbose) print(result)
  invisible(result)
}