# Revised statistical inference functions for empirical Markov chains.

#' Test the first-order Markov property of an empirical sequence
#'
#' Tests the null hypothesis that the conditional distribution of the next
#' state depends only on the current state, against a second-order Markov
#' alternative. Equivalently, it tests conditional independence of the
#' previous and next states given the current state.
#'
#' @rdname statisticalTests
#' @family statisticalTests
#' @param sequence An empirical sequence of states.
#' @param method Test statistic: `"G"`, `"Pearson"`, or `"simulation"`.
#' @param B Number of Monte Carlo replicates for `method = "simulation"`.
#' @param seed Optional random seed.
#' @param verbose Should results be printed?
#' @return A list with statistic, degrees of freedom, p-value, observed and
#' expected triple-transition counts, and the fitted first-order transition
#' matrix. For simulation, simulated statistics and `B` are also returned.
#' @references Anderson, T. W. and Goodman, L. A. (1957). Statistical
#' inference about Markov chains. *The Annals of Mathematical Statistics*,
#' 28(1), 89--110.
#' @export
verifyMarkovProperty <- function(sequence,
                                 method = c("G", "Pearson", "simulation"),
                                 B = 9999, seed = NULL, verbose = TRUE) {
  method <- match.arg(method)
  if (length(sequence) < 4L) stop("sequence must contain at least four observations.")
  if (anyNA(sequence)) stop("sequence must not contain missing values.")
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B))
    stop("B must be a positive integer.")

  states <- unique(sequence)
  r <- length(states)
  n <- length(sequence)
  idx <- match(sequence, states)

  Nijk <- array(0L, dim = c(r, r, r), dimnames = list(states, states, states))
  linear <- idx[-c(1L, 2L)] +
    (idx[2:(n - 1L)] - 1L) * r +
    (idx[1:(n - 2L)] - 1L) * r * r
  Nijk[] <- tabulate(linear, nbins = r^3L)

  Nij <- apply(Nijk, c(1L, 2L), sum)
  Njk <- apply(Nijk, c(2L, 3L), sum)
  Nj <- rowSums(Njk)
  expected <- array(0, dim = c(r, r, r), dimnames = dimnames(Nijk))
  for (j in seq_len(r)) if (Nj[j] > 0L)
    expected[, j, ] <- outer(Nij[, j], Njk[j, ]) / Nj[j]

  statistic_fun <- function(observed, expected, method) {
    positive_expected <- expected > 0
    if (method == "Pearson")
      return(sum((observed[positive_expected] - expected[positive_expected])^2 /
                 expected[positive_expected]))
    positive <- positive_expected & observed > 0
    2 * sum(observed[positive] * log(observed[positive] / expected[positive]))
  }

  statistic <- statistic_fun(Nijk, expected,
                             if (method == "Pearson") "Pearson" else "G")
  dof <- 0L
  for (j in seq_len(r)) {
    previous <- sum(Nij[, j] > 0)
    future <- sum(Njk[j, ] > 0)
    if (previous > 0L && future > 0L)
      dof <- dof + (previous - 1L) * (future - 1L)
  }

  P <- matrix(0, r, r, dimnames = list(states, states))
  informative <- Nj > 0
  P[informative, ] <- sweep(Njk[informative, , drop = FALSE], 1L,
                            Nj[informative], "/")

  if (method != "simulation") {
    p.value <- if (dof > 0) pchisq(statistic, dof, lower.tail = FALSE) else NA_real_
    result <- list(statistic = statistic, dof = dof, p.value = p.value,
                   method = method, observed = Nijk, expected = expected,
                   transitionMatrix = P)
  } else {
    if (!all(informative)) stop("method = 'simulation' requires at least one observed outgoing transition from every state appearing in the sequence.")
    if (!is.null(seed)) set.seed(seed)
    simulated <- numeric(B)
    for (b in seq_len(B)) {
      x <- integer(n); x[1L] <- idx[1L]
      for (t in seq_len(n - 1L)) x[t + 1L] <- sample.int(r, 1L, prob = P[x[t], ])
      a <- x[-c(1L, 2L)] + (x[2:(n - 1L)] - 1L) * r +
        (x[1:(n - 2L)] - 1L) * r * r
      Ns <- array(tabulate(a, nbins = r^3L), dim = c(r, r, r))
      Nijs <- apply(Ns, c(1L, 2L), sum)
      Njks <- apply(Ns, c(2L, 3L), sum)
      Njs <- rowSums(Njks)
      Es <- array(0, dim = c(r, r, r))
      for (j in seq_len(r)) if (Njs[j] > 0L)
        Es[, j, ] <- outer(Nijs[, j], Njks[j, ]) / Njs[j]
      simulated[b] <- statistic_fun(Ns, Es,
                                     if (method == "Pearson") "Pearson" else "G")
    }
    p.value <- (1 + sum(simulated >= statistic)) / (B + 1)
    result <- list(statistic = statistic, dof = dof, p.value = p.value,
                   method = method, observed = Nijk, expected = expected,
                   transitionMatrix = P, simulations = simulated, B = B,
                   initial.state = states[idx[1L]])
  }
  if (verbose) cat("Testing compatibility with a first-order Markov model\n",
                   "Method:", method, "\nStatistic:", statistic,
                   "\nDegrees of freedom:", dof, "\np-value:", p.value, "\n")
  invisible(result)
}

#' Assess the order of an empirical Markov chain
#'
#' Tests the first-order Markov property state by state by comparing the
#' distribution of the state two steps ahead across possible predecessor
#' states. The statistic is the sum of Pearson chi-squared statistics.
#'
#' @rdname statisticalTests
#' @family statisticalTests
#' @param sequence An empirical sequence of states.
#' @param verbose Should test results be printed?
#' @return A list containing the statistic, degrees of freedom and p-value.
#' @export
assessOrder <- function(sequence, verbose = TRUE) {
  if (length(sequence) < 4L) stop("sequence must contain at least four observations.")
  if (anyNA(sequence)) stop("sequence must not contain missing values.")

  states <- unique(sequence)
  k <- length(states)
  n <- length(sequence)
  statistic <- 0

  for (present in states) {
    mat <- matrix(0, nrow = k, ncol = k,
                  dimnames = list(states, states))
    for (i in seq_len(n - 2L)) {
      if (identical(sequence[i + 1L], present)) {
        past <- as.character(sequence[i])
        future <- as.character(sequence[i + 2L])
        mat[past, future] <- mat[past, future] + 1
      }
    }

    row_totals <- rowSums(mat)
    active <- row_totals > 0
    if (sum(active) > 1L && sum(colSums(mat) > 0) > 1L) {
      expected <- outer(row_totals[active], colSums(mat)) /
        sum(row_totals[active])
      observed <- mat[active, , drop = FALSE]
      positive <- expected > 0
      statistic <- statistic +
        sum((observed[positive] - expected[positive])^2 / expected[positive])
    }
  }

  dof <- k * (k - 1L)^2
  p.value <- if (dof > 0L) pchisq(statistic, dof, lower.tail = FALSE) else NA_real_
  result <- list(statistic = statistic, dof = dof, p.value = p.value)

  if (verbose) {
    cat("The assessOrder test statistic is:", statistic, "\n",
        "The Chi-Square d.f. are:", dof, "\n",
        "The p-value is:", p.value, "\n")
  }
  invisible(result)
}

#' Test compatibility of an empirical transition matrix with a theoretical one
#'
#' @rdname statisticalTests
#' @family statisticalTests
#' @param data An empirical sequence or a matrix of transition counts.
#' @param object A `markovchain` object specifying theoretical probabilities.
#' @param method Test statistic: `"G"`, `"Pearson"`, or `"simulation"`.
#' @param B Number of Monte Carlo replicates for `method = "simulation"`.
#' @param seed Optional random seed.
#' @param verbose Should results be printed?
#' @return A list containing statistic, degrees of freedom, p-value, method,
#' observed counts, expected counts, and theoretical transition matrix.
#' @references Kullback, S., Kupperman, M. and Ku, H. H. (1962). Tests for
#' Contingency Tables and Markov Chains. *Technometrics*, 4(4), 573--608.
#' @export
verifyEmpiricalToTheoretical <- function(data, object,
                                         method = c("G", "Pearson", "simulation"),
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
  observed <- observed[states, states, drop = FALSE]
  row_totals <- rowSums(observed)
  expected <- sweep(P, 1L, row_totals, "*")
  informative <- row_totals > 0
  if (any(observed > 0 & P == 0)) {
    result <- list(statistic = Inf, dof = NA_integer_, p.value = 0, method = method,
                   observed = observed, expected = expected, transitionMatrix = P)
    if (verbose) cat("At least one observed transition has zero theoretical probability; the null hypothesis is rejected.\n")
    return(invisible(result))
  }
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
    if (!is.null(seed)) set.seed(seed)
    simulated <- numeric(B)
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
  if (verbose) cat("Testing compatibility of empirical and theoretical transition probabilities\n",
                   "Method:", method, "\nStatistic:", statistic,
                   "\nDegrees of freedom:", dof, "\np-value:", p.value, "\n")
  invisible(result)
}

#' Test homogeneity of transition probabilities across empirical sequences
#'
#' @rdname statisticalTests
#' @family statisticalTests
#' @param inputList A list whose elements are empirical sequences or matrices.
#' @param method Test statistic: `"G"`, `"Pearson"`, or `"simulation"`.
#' @param B Number of Monte Carlo replicates for `method = "simulation"`.
#' @param seed Optional random seed.
#' @param verbose Should results be printed?
#' @return A list containing statistic, degrees of freedom, p-value, method,
#' pooled transition counts, and individual transition counts.
#' @references Kullback, S., Kupperman, M. and Ku, H. H. (1962). Tests for
#' Contingency Tables and Markov Chains. *Technometrics*, 4(4), 573--608.
#' @export
verifyHomogeneity <- function(inputList, method = c("G", "Pearson", "simulation"),
                              B = 9999, seed = NULL, verbose = TRUE) {
  method <- match.arg(method)
  if (!is.list(inputList) || length(inputList) < 2L) stop("inputList must contain at least two sequences or matrices.")
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B)) stop("B must be a positive integer.")
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
    out <- matrix(0, length(states), length(states), dimnames = list(states, states))
    out[rownames(m), colnames(m)] <- m; out
  })
  pooled <- Reduce("+", mats)
  hom_stat <- function(mats0) {
    pooled0 <- Reduce("+", mats0); stat <- 0
    for (i in seq_along(states)) {
      totals <- vapply(mats0, function(m) sum(m[i, ]), numeric(1)); active <- totals > 0
      if (sum(active) < 2L) next
      destinations <- pooled0[i, ] > 0; k <- sum(destinations)
      if (k < 2L) next
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
    p.value <- if (dof > 0) pchisq(statistic, dof, lower.tail = FALSE) else NA_real_; simulated <- NULL
  } else {
    if (!is.null(seed)) set.seed(seed)
    simulated <- numeric(B)
    for (b in seq_len(B)) {
      sim_mats <- lapply(mats, function(m) matrix(0, nrow(m), ncol(m), dimnames = dimnames(m)))
      pooled0 <- pooled
      for (i in seq_along(states)) {
        destinations <- pooled0[i, ] > 0; if (!any(destinations)) next
        p <- pooled0[i, destinations] / sum(pooled0[i, destinations])
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
  if (verbose) cat("Testing homogeneity of transition probabilities across sequences\n",
                   "Method:", method, "\nStatistic:", statistic,
                   "\nDegrees of freedom:", dof, "\np-value:", p.value, "\n")
  invisible(result)
}
