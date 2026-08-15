#' Test compatibility of an empirical transition matrix with a theoretical one
#'
#' Tests whether observed transition counts are compatible with the transition
#' probabilities of a theoretical `markovchain` object. The null hypothesis is
#' that, conditional on the current state, the next state follows the
#' probabilities specified by `object`.
#'
#' The default likelihood-ratio statistic is
#' `2 * sum(n_ij log(n_ij / (n_i+ p_ij)))`. A Pearson statistic is also
#' available. Rows with no observed transitions are not informative.
#'
#' @rdname statisticalTests
#' @family statisticalTests
#' @param data An empirical sequence or a matrix of transition counts. A
#'   sequence is converted to transition counts with `createSequenceMatrix()`.
#' @param object A `markovchain` object specifying the theoretical transition
#'   probabilities.
#' @param method Test statistic: `"G"` or `"Pearson"`.
#' @param verbose Should results be printed?
#' @return A list containing `statistic`, `dof`, `p.value`, `method`, observed
#'   counts, expected counts, and the theoretical transition matrix.
#' @references Kullback, S., Kupperman, M. and Ku, H. H. (1962). Tests for
#'   Contingency Tables and Markov Chains. *Technometrics*, 4(4), 573--608.
#' @export
verifyEmpiricalToTheoretical <- function(data, object,
                                         method = c("G", "Pearson"),
                                         verbose = TRUE) {
  method <- match.arg(method)
  if (missing(data) || missing(object))
    stop("Both data and object must be supplied.")
  if (!methods::is(object, "markovchain"))
    stop("object must be a markovchain object.")

  states <- names(object)
  P <- object@transitionMatrix
  if (is.null(rownames(P)) || is.null(colnames(P)))
    stop("The theoretical transition matrix must have state names.")

  if (is.matrix(data)) {
    observed <- data
    if (is.null(rownames(observed)) || is.null(colnames(observed)))
      stop("A transition-count matrix must have row and column names.")
  } else if (is.atomic(data)) {
    observed <- createSequenceMatrix(stringchar = data, possibleStates = states)
  } else {
    stop("data must be an empirical sequence or a transition-count matrix.")
  }

  if (anyNA(observed) || any(observed < 0) || any(observed != floor(observed)))
    stop("Observed transition counts must be non-negative integers.")
  if (!setequal(rownames(observed), states) || !setequal(colnames(observed), states))
    stop("Empirical and theoretical transition matrices must have the same state support.")

  observed <- observed[states, states, drop = FALSE]
  expected <- sweep(P, 1L, rowSums(observed), "*")
  informative <- rowSums(observed) > 0

  if (any(observed[informative, , drop = FALSE] > 0 & P[informative, , drop = FALSE] == 0)) {
    result <- list(statistic = Inf, dof = NA_integer_, p.value = 0,
                   method = method, observed = observed,
                   expected = expected, transitionMatrix = P)
    if (verbose)
      cat("At least one observed transition has zero theoretical probability; the null hypothesis is rejected.\n")
    return(invisible(result))
  }

  positive <- informative & P > 0
  if (method == "G") {
    statistic <- 2 * sum(observed[positive] *
                          log(observed[positive] / expected[positive]))
  } else {
    statistic <- sum((observed[positive] - expected[positive])^2 /
                     expected[positive])
  }

  dof <- sum(vapply(which(informative), function(i)
    sum(P[i, ] > 0) - 1L, integer(1)))
  p.value <- if (dof > 0) pchisq(statistic, df = dof, lower.tail = FALSE) else NA_real_

  result <- list(statistic = unname(statistic), dof = dof,
                 p.value = p.value, method = method,
                 observed = observed, expected = expected,
                 transitionMatrix = P)
  if (verbose) {
    cat("Testing compatibility of empirical and theoretical transition probabilities\n")
    cat("Method:", method, "\n")
    cat("Statistic:", statistic, "\n")
    cat("Degrees of freedom:", dof, "\n")
    cat("p-value:", p.value, "\n")
  }
  invisible(result)
}

#' Test homogeneity of transition probabilities across empirical sequences
#'
#' Tests whether several empirical transition matrices can be regarded as
#' realizations of the same time-homogeneous Markov chain. The null hypothesis
#' is that, for every origin state, the conditional distribution of the next
#' state is the same in all input sequences.
#'
#' @rdname statisticalTests
#' @family statisticalTests
#' @param inputList A list whose elements are empirical sequences or matrices
#'   of transition counts.
#' @param method Test statistic: `"G"` or `"Pearson"`.
#' @param verbose Should results be printed?
#' @return A list containing `statistic`, `dof`, `p.value`, `method`, the
#'   pooled transition counts, and the individual transition-count matrices.
#' @references Kullback, S., Kupperman, M. and Ku, H. H. (1962). Tests for
#'   Contingency Tables and Markov Chains. *Technometrics*, 4(4), 573--608.
#' @export
verifyHomogeneity <- function(inputList,
                              method = c("G", "Pearson"),
                              verbose = TRUE) {
  method <- match.arg(method)
  if (!is.list(inputList) || length(inputList) < 2L)
    stop("inputList must be a list containing at least two sequences or matrices.")

  mats <- lapply(inputList, function(x) {
    if (is.matrix(x)) {
      if (is.null(rownames(x)) || is.null(colnames(x)))
        stop("Transition-count matrices must have row and column names.")
      m <- x
    } else if (is.atomic(x)) {
      m <- createSequenceMatrix(stringchar = x)
    } else {
      stop("Each element of inputList must be a sequence or transition-count matrix.")
    }
    if (anyNA(m) || any(m < 0) || any(m != floor(m)))
      stop("Transition counts must be non-negative integers.")
    m
  })

  states <- sort(unique(unlist(lapply(mats, rownames))))
  mats <- lapply(mats, function(m) {
    out <- matrix(0, nrow = length(states), ncol = length(states),
                  dimnames = list(states, states))
    out[rownames(m), colnames(m)] <- m
    out
  })

  pooled <- Reduce("+", mats)
  statistic <- 0
  dof <- 0L

  for (i in seq_along(states)) {
    totals <- vapply(mats, function(m) sum(m[i, ]), numeric(1))
    active <- totals > 0
    if (sum(active) < 2L)
      next

    pooledRow <- pooled[i, ]
    destinations <- pooledRow > 0
    k <- sum(destinations)
    if (k < 2L)
      next

    p <- pooledRow[destinations] / sum(pooledRow)
    for (s in which(active)) {
      obs <- mats[[s]][i, destinations]
      exp <- totals[s] * p
      if (method == "G") {
        positive <- obs > 0
        statistic <- statistic + 2 * sum(obs[positive] *
          log(obs[positive] / exp[positive]))
      } else {
        statistic <- statistic + sum((obs - exp)^2 / exp)
      }
    }
    dof <- dof + (sum(active) - 1L) * (k - 1L)
  }

  p.value <- if (dof > 0) pchisq(statistic, df = dof, lower.tail = FALSE) else NA_real_
  result <- list(statistic = unname(statistic), dof = dof,
                 p.value = p.value, method = method,
                 pooled = pooled, transitionCounts = mats)
  if (verbose) {
    cat("Testing homogeneity of transition probabilities across sequences\n")
    cat("Method:", method, "\n")
    cat("Statistic:", statistic, "\n")
    cat("Degrees of freedom:", dof, "\n")
    cat("p-value:", p.value, "\n")
  }
  invisible(result)
}
