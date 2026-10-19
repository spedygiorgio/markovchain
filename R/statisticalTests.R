# Revised statistical inference functions for empirical Markov chains.

#' Standardize statistical inference results to R's htest convention.
#'
#' @keywords internal
# Converts a plain list (statistic, dof, p.value, ...) into the standard R 'htest' object, so that
# print() shows the usual test layout (statistic, degrees of freedom, p-value). `method` selects the
# labels shown to the user.
.asHtest <- function(result, method, data.name) {
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

#' Statistical tests for Markov chains
#'
#' @description Tests of hypotheses on observed sequences of states, based on
#'   transition counts (Anderson and Goodman 1957; Kullback, Kupperman and Ku
#'   1962):
#'   * `verifyMarkovProperty()`: the first-order Markov property, i.e. that
#'     the next state is independent of the previous one given the present one
#'     (against a second-order chain);
#'   * `assessOrder()`: the same hypothesis with the Pearson statistic, also
#'     for a list of independent sequences;
#'   * `verifyEmpiricalToTheoretical()`: that the sequence (or a matrix of
#'     transition counts) comes from a given chain;
#'   * `verifyHomogeneity()`: that several sequences (or count matrices) share
#'     the same transition matrix;
#'   * `assessStationarity()`: that the transition probabilities do not change
#'     over time, comparing `nblocks` consecutive blocks of the sequence.
#'
#' @details Each test compares observed transition counts with the counts
#'   expected under the null hypothesis, by the likelihood-ratio statistic
#'   \eqn{G^2 = 2 \sum O \log(O / E)} (`method = "G"`) or by Pearson's
#'   \eqn{X^2 = \sum (O - E)^2 / E} (`method = "Pearson"`); cells with no
#'   expected count are left out, and the degrees of freedom count only the
#'   transitions actually observed. The p-value is from the chi-squared
#'   distribution, or, with `method = "simulation"`, from `B` samples
#'   simulated under the null hypothesis (the \eqn{G^2} statistic, p-value
#'   \eqn{(1 + \#\{simulated \ge observed\}) / (B + 1)}), which is more
#'   reliable when many counts are small. Sequences of numbers or factors are
#'   compared as character strings.
#'
#' @return An object of class `"htest"`, printed if `verbose = TRUE` and
#'   returned invisibly, with the usual components `statistic`, `parameter`
#'   (degrees of freedom), `p.value`, `method` and `data.name`, and also
#'   `dof`. Depending on the test it also contains the `observed` and
#'   `expected` counts and the estimated `transitionMatrix`
#'   (`verifyMarkovProperty()`, `verifyEmpiricalToTheoretical()`), the
#'   `pooled` counts and the `transitionCounts` of each element
#'   (`verifyHomogeneity()`), and, with `method = "simulation"`, the
#'   simulated statistics (`simulations`) and `B`.
#'
#' @references
#' Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about
#' Markov chains. The Annals of Mathematical Statistics 28(1), 89-110.
#'
#' Kullback, S., Kupperman, M. and Ku, H. H. (1962). Tests for contingency
#' tables and Markov chains. Technometrics 4(4), 573-608.
#'
#' @seealso [assessIndependence()] (independence of consecutive states),
#'   [selectOrder()] (order by information criteria), [compareChains()].
#' @examples
#' data(rain, package = "markovchain")
#' # first-order Markov property of the Alofi rainfall sequence
#' verifyMarkovProperty(rain$rain)
#' assessOrder(rain$rain)
#'
#' # compatibility with a chain: the fitted one, and a different one
#' fit <- markovchainFit(rain$rain)$estimate
#' verifyEmpiricalToTheoretical(rain$rain, fit)
#' s <- c("0", "1-5", "6+")
#' uniform <- new("markovchain", states = s,
#'                transitionMatrix = matrix(1 / 3, 3, 3, dimnames = list(s, s)))
#' verifyEmpiricalToTheoretical(rain$rain, uniform)
#'
#' # the first and second half of the sequence have the same transitions
#' verifyHomogeneity(list(rain$rain[1:548], rain$rain[549:1096]))
#' # time-homogeneity over four blocks
#' assessStationarity(rain$rain, nblocks = 4)
#'
#' # a Monte Carlo p-value
#' verifyMarkovProperty(rain$rain, method = "simulation", B = 499, seed = 1)
#' @name statisticalTests
NULL

#' @rdname statisticalTests
#' @param sequence An observed sequence of states (a vector); for
#'   `assessOrder()` also a list of independent sequences, whose triples of
#'   consecutive states are formed within each sequence and pooled.
#' @param method The statistic: `"G"` (likelihood ratio, the default),
#'   `"Pearson"`, or `"simulation"` (likelihood ratio with a Monte Carlo
#'   p-value).
#' @param B Number of Monte Carlo samples for `method = "simulation"`.
#' @param seed Optional random seed for `method = "simulation"`.
#' @param verbose Logical: print the result?
#' @export
# Test of the Markov property (Anderson-Goodman): is the next state independent of the PAST given the
# PRESENT? For every present state j it builds the table of (past i, future k) counts and compares it
# with the counts expected under independence.
verifyMarkovProperty <- function(sequence, method = c("G", "Pearson", "simulation"),
                                 B = 9999, seed = NULL, verbose = TRUE) {
  method <- match.arg(method)
  if (length(sequence) < 4L) stop("sequence must contain at least four observations.")
  if (anyNA(sequence)) stop("sequence must not contain missing values.")
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B))
    stop("B must be a positive integer.")
  # Step 1: count triples (i, j, k) of consecutive observations in a 3D array Nijk;
  # `linear` converts each triple into one index of the array so tabulate() can count them quickly.
  states <- unique(sequence); r <- length(states); n <- length(sequence)
  idx <- match(sequence, states)
  Nijk <- array(0L, dim = c(r, r, r), dimnames = list(states, states, states))
  linear <- idx[-c(1L, 2L)] + (idx[2:(n - 1L)] - 1L) * r +
    (idx[1:(n - 2L)] - 1L) * r * r
  Nijk[] <- tabulate(linear, nbins = r^3L)
  # Marginals: Nij (pairs past-present), Njk (pairs present-future), Nj (times in state j with a future).
  Nij <- apply(Nijk, c(1L, 2L), sum); Njk <- apply(Nijk, c(2L, 3L), sum)
  Nj <- rowSums(Njk)
  # Step 2: expected counts under the Markov property: N_ij * N_jk / N_j.
  expected <- array(0, dim = c(r, r, r), dimnames = dimnames(Nijk))
  for (j in seq_len(r)) if (Nj[j] > 0L)
    expected[, j, ] <- outer(Nij[, j], Njk[j, ]) / Nj[j]
  # Step 3: test statistic, Pearson (sum (O-E)^2/E) or G (2 sum O log(O/E)); cells with zero
  # expectation or zero observation are skipped to avoid dividing by zero or taking log(0).
  statistic_fun <- function(observed, expected, method) {
    positive_expected <- expected > 0
    if (method == "Pearson")
      return(sum((observed[positive_expected] - expected[positive_expected])^2 /
                 expected[positive_expected]))
    positive <- positive_expected & observed > 0
    2 * sum(observed[positive] * log(observed[positive] / expected[positive]))
  }
  statistic <- statistic_fun(Nijk, expected, if (method == "Pearson") "Pearson" else "G")
  # Degrees of freedom: for each present state, (#observed pasts - 1) * (#observed futures - 1).
  dof <- 0L
  for (j in seq_len(r)) {
    previous <- sum(Nij[, j] > 0); future <- sum(Njk[j, ] > 0)
    if (previous > 0L && future > 0L) dof <- dof + (previous - 1L) * (future - 1L)
  }
  # Estimated transition matrix: P[j,k] = Njk / Nj.
  P <- matrix(0, r, r, dimnames = list(states, states)); informative <- Nj > 0
  P[informative, ] <- sweep(Njk[informative, , drop = FALSE], 1L, Nj[informative], "/")
  # Asymptotic p-value from the chi-squared distribution ...
  if (method != "simulation") {
    p.value <- if (dof > 0) pchisq(statistic, dof, lower.tail = FALSE) else NA_real_
    result <- list(statistic = statistic, dof = dof, p.value = p.value, method = method,
                   observed = Nijk, expected = expected, transitionMatrix = P)
  # ... or a Monte Carlo p-value: simulate B sequences from the estimated chain, recompute the statistic
  # on each, and report the share of simulated statistics at least as large as the observed one.
  } else {
    if (!all(informative)) stop("method = 'simulation' requires at least one observed outgoing transition from every state appearing in the sequence.")
    if (!is.null(seed)) set.seed(seed)
    # The B sequences are simulated in C++ (src/markovPropertyTest.cpp), each
    # step drawn exactly as sample.int(r, 1, prob = P[x, ]) would draw it, and
    # the G statistic of each is computed with the same operations as
    # statistic_fun(): same random numbers as an R loop, and the same values up
    # to the last bits (bit for bit on x86-64). The observed statistic is
    # recomputed by the same C++ function, so that a simulated sequence with
    # the observed triple counts gives exactly the observed value on every
    # platform and counts as a tie (>=) in the p-value.
    storage.mode(P) <- "double"
    simulated <- .markovPropertySimulationRcpp(P, idx[1L], n, as.integer(B), FALSE)
    statistic <- .markovPropertyStatisticRcpp(idx, r, FALSE)
    p.value <- (1 + sum(simulated >= statistic)) / (B + 1)
    result <- list(statistic = statistic, dof = dof, p.value = p.value, method = method,
                   observed = Nijk, expected = expected, transitionMatrix = P,
                   simulations = simulated, B = B, initial.state = states[idx[1L]])
  }
  # Wrap the result into an 'htest' object.
  result <- .asHtest(result, method, deparse(substitute(sequence)))
  if (verbose) print(result)
  invisible(result)
}

#' @rdname statisticalTests
#' @export
# Tests the order of the chain: is the future independent of the past given the present? (the same idea
# as verifyMarkovProperty, but with the Pearson statistic pooled over present states.)
assessOrder <- function(sequence, verbose = TRUE) {
  data.name <- deparse(substitute(sequence))
  # A list means several independent sequences: triples are taken inside each sequence and pooled.
  # a list is a set of independent sequences: the triples (past, present,
  # future) are taken within each sequence and pooled
  if (is.list(sequence)) {
    if (length(sequence) == 0L || !all(vapply(sequence, is.atomic, logical(1)))) {
      stop("sequence must be a vector of states or a list of such vectors.")
    }
    if (anyNA(unlist(sequence, use.names = FALSE))) stop("sequence must not contain missing values.")
    seqs <- lapply(sequence, as.character)
    if (sum(pmax(lengths(seqs) - 2L, 0L)) < 2L) {
      stop("sequence must contain at least four observations.")
    }
  } else {
    if (length(sequence) < 4L) stop("sequence must contain at least four observations.")
    if (anyNA(sequence)) stop("sequence must not contain missing values.")
    # factors (and numbers) are compared as character, like in verifyMarkovProperty
    seqs <- list(as.character(sequence))
  }
  seqs <- seqs[lengths(seqs) >= 3L]
  states <- unique(unlist(seqs, use.names = FALSE))
  past <- factor(unlist(lapply(seqs, function(s) s[seq_len(length(s) - 2L)])), levels = states)
  present_state <- unlist(lapply(seqs, function(s) s[2:(length(s) - 1L)]))
  future <- factor(unlist(lapply(seqs, function(s) s[3:length(s)])), levels = states)
  # For each present state, cross-tabulate the past against the future and test independence.
  statistic <- 0; dof <- 0L
  for (present in states) {
    selected <- present_state == present
    mat <- unclass(table(past[selected], future[selected]))
    row_totals <- rowSums(mat); active <- row_totals > 0
    active_cols <- sum(colSums(mat) > 0)
    if (sum(active) > 1L && active_cols > 1L) {
      expected <- outer(row_totals[active], colSums(mat)) / sum(row_totals[active])
      observed <- mat[active, , drop = FALSE]; positive <- expected > 0
      statistic <- statistic + sum((observed[positive] - expected[positive])^2 / expected[positive])
      dof <- dof + (sum(active) - 1L) * (active_cols - 1L)
    }
  }
  result <- list(statistic = statistic, dof = dof,
                 p.value = if (dof > 0L) pchisq(statistic, dof, lower.tail = FALSE) else NA_real_)
  result <- .asHtest(result, "Pearson", data.name)
  result$method <- "Pearson's Chi-squared test for Markov order"
  if (verbose) print(result)
  invisible(result)
}

#' @rdname statisticalTests
#' @param data An observed sequence of states, or a matrix of transition
#'   counts with the states as dimnames.
#' @param object A [`markovchain`][markovchain-class] object, the chain of the
#'   null hypothesis.
#' @export
# Tests whether an observed sequence (or a matrix of transition counts) is compatible with a given
# theoretical chain `object`. Expected counts are rowTotals[i] * P[i, j].
verifyEmpiricalToTheoretical <- function(data, object, method = c("G", "Pearson", "simulation"),
                                         B = 9999, seed = NULL, verbose = TRUE) {
  method <- match.arg(method)
  if (missing(data) || missing(object)) stop("Both data and object must be supplied.")
  if (!methods::is(object, "markovchain")) stop("object must be a markovchain object.")
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B)) stop("B must be a positive integer.")
  states <- names(object); P <- object@transitionMatrix
  if (is.matrix(data)) {
    observed <- data
    if (is.null(rownames(observed)) || is.null(colnames(observed))) stop("A transition-count matrix must have row and column names.")
  } else if (is.atomic(data)) observed <- createSequenceMatrix(stringchar = data, possibleStates = states)
  else stop("data must be an empirical sequence or a transition-count matrix.")
  if (anyNA(observed) || any(observed < 0) || any(observed != floor(observed))) stop("Observed transition counts must be non-negative integers.")
  if (!setequal(rownames(observed), states) || !setequal(colnames(observed), states)) stop("Empirical and theoretical transition matrices must have the same state support.")
  observed <- observed[states, states, drop = FALSE]; row_totals <- rowSums(observed)
  expected <- sweep(P, 1L, row_totals, "*"); informative <- row_totals > 0
  # Impossible transition observed (probability 0 in the theory): the data cannot come from the chain.
  if (any(observed > 0 & P == 0)) {
    result <- list(statistic = Inf, dof = NA_integer_, p.value = 0, method = method,
                   observed = observed, expected = expected, transitionMatrix = P)
  } else {
    # Test statistic (Pearson or G).
    statistic <- if (method == "Pearson") {
      pos <- outer(informative, rep(TRUE, length(states))) & P > 0
      sum((observed[pos] - expected[pos])^2 / expected[pos])
    } else {
      pos <- P > 0 & observed > 0
      2 * sum(observed[pos] * log(observed[pos] / expected[pos]))
    }
    # Degrees of freedom: for each row, the number of possible transitions minus one.
    dof <- sum(vapply(which(informative), function(i) sum(P[i, ] > 0) - 1L, integer(1)))
    if (method != "simulation") {
      p.value <- if (dof > 0) pchisq(statistic, dof, lower.tail = FALSE) else NA_real_
      result <- list(statistic = statistic, dof = dof, p.value = p.value, method = method,
                     observed = observed, expected = expected, transitionMatrix = P)
    # Monte Carlo version: simulate the transition counts from a multinomial distribution row by row.
    } else {
      if (!is.null(seed)) set.seed(seed); simulated <- numeric(B)
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

#' @rdname statisticalTests
#' @param inputList A list of at least two observed sequences or matrices of
#'   transition counts (with the states as dimnames).
#' @export
# Tests whether several sequences (or count matrices) share the same transition matrix:
# for each state the rows of the different samples are compared with the pooled row.
verifyHomogeneity <- function(inputList, method = c("G", "Pearson", "simulation"),
                              B = 9999, seed = NULL, verbose = TRUE) {
  method <- match.arg(method)
  if (!is.list(inputList) || length(inputList) < 2L) stop("inputList must contain at least two sequences or matrices.")
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B)) stop("B must be a positive integer.")
  # Step 1: obtain a matrix of transition counts from each element of the list.
  mats <- lapply(inputList, function(x) {
    if (is.matrix(x)) {
      if (is.null(rownames(x)) || is.null(colnames(x))) stop("Transition-count matrices must have row and column names.")
      m <- x
    } else if (is.atomic(x)) m <- createSequenceMatrix(stringchar = x)
    else stop("Each element must be a sequence or transition-count matrix.")
    if (anyNA(m) || any(m < 0) || any(m != floor(m))) stop("Transition counts must be non-negative integers.")
    m
  })
  # Step 2: put all matrices on the same set of states (missing states get zero counts).
  states <- sort(unique(unlist(lapply(mats, rownames))))
  mats <- lapply(mats, function(m) {
    out <- matrix(0, length(states), length(states), dimnames = list(states, states)); out[rownames(m), colnames(m)] <- m; out
  })
  # Step 3: pooled counts = the sum of all samples.
  pooled <- Reduce("+", mats)
  # hom_stat(): statistic that compares, state by state, each sample's counts with the counts expected
  # if that sample used the pooled probabilities.
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
  # Degrees of freedom: (#samples with data - 1) * (#possible destinations - 1) for each state.
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

#' @rdname statisticalTests
#' @param nblocks Number of consecutive blocks of (nearly) equal length into
#'   which the sequence is split, at least two.
#' @param structural.zeros Optional logical matrix, with one row and column
#'   per state (in the order of `unique(sequence)`, or with the states as
#'   dimnames), marking the transitions that are impossible; they are left
#'   out of the comparison.
#' @export
# Tests time-homogeneity: split the sequence into `nblocks` blocks and check, for every state, whether
# the distribution of the next state is the same in all blocks (chi-squared test of independence on the
# block x next-state table). `structural.zeros` marks transitions known to be impossible.
assessStationarity <- function(sequence, nblocks, structural.zeros = NULL, verbose = TRUE) {
  if (length(sequence) < 2L) stop("sequence must contain at least two observations")
  if (length(nblocks) != 1L || !is.numeric(nblocks) || !is.finite(nblocks) ||
      nblocks < 2 || nblocks != as.integer(nblocks))
    stop("nblocks must be a positive integer greater than one")
  # Step 1: validate the arguments and the optional matrix of structural zeros.
  nblocks <- as.integer(nblocks); states <- unique(sequence); nstates <- length(states)
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
  # Which transitions occur at least once in the data.
  observed <- matrix(FALSE, nstates, nstates, dimnames = list(as.character(states), as.character(states)))
  for (j in seq_len(length(sequence) - 1L)) {
    from <- match(sequence[j], states); to <- match(sequence[j + 1L], states); observed[from, to] <- TRUE
  }
  if (any(observed & structural.zeros)) stop("sequence contains a transition marked as a structural zero")
  # Step 2: for every state i, count the transitions out of i in each block (the block is derived from
  # the position j in the sequence), drop empty blocks and structurally impossible columns,
  # and accumulate Pearson's statistic and its degrees of freedom.
  blocksize <- length(sequence) / nblocks; TStat <- 0; df <- 0
  for (i in seq_len(nstates)) {
    counts <- matrix(0, nrow = nblocks, ncol = nstates, dimnames = list(seq_len(nblocks), as.character(states)))
    for (j in seq_len(length(sequence) - 1L)) if (sequence[j] == states[i]) {
      block <- min(nblocks, ceiling(j / blocksize)); future <- match(sequence[j + 1L], states); counts[block, future] <- counts[block, future] + 1
    }
    counts <- counts[rowSums(counts) > 0, , drop = FALSE]; counts <- counts[, !structural.zeros[i, ], drop = FALSE]
    if (nrow(counts) < 2L || ncol(counts) < 2L) next
    expected <- outer(rowSums(counts), colSums(counts) / sum(counts)); positive.expected <- expected > 0
    TStat <- TStat + sum((counts[positive.expected] - expected[positive.expected])^2 / expected[positive.expected])
    df <- df + (nrow(counts) - 1L) * (ncol(counts) - 1L)
  }
  # Step 3: p-value from the chi-squared distribution.
  pvalue <- if (df > 0) pchisq(TStat, df = df, lower.tail = FALSE) else NA_real_
  result <- list(statistic = TStat, dof = df, p.value = pvalue)
  result <- .asHtest(result, "Pearson", deparse(substitute(sequence)))
  result$method <- "Pearson's Chi-squared test for time-homogeneity"
  if (verbose) print(result)
  invisible(result)
}
