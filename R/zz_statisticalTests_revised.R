#' Test compatibility of an empirical transition matrix with a theoretical one
#'
#' Tests whether observed transition counts are compatible with the transition
#' probabilities of a theoretical `markovchain` object.
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
#'   observed counts, expected counts, and theoretical transition matrix.
#'   For simulation, the simulated statistics and number of replicates are
#'   also returned.
#' @references Kullback, S., Kupperman, M. and Ku, H. H. (1962). Tests for
#'   Contingency Tables and Markov Chains. *Technometrics*, 4(4), 573--608.
#' @export
verifyEmpiricalToTheoretical <- function(data, object,
                                         method = c("G", "Pearson", "simulation"),
                                         B = 9999, seed = NULL,
                                         verbose = TRUE) {
  method <- match.arg(method)
  if (missing(data) || missing(object)) stop("Both data and object must be supplied.")
  if (!methods::is(object, "markovchain")) stop("object must be a markovchain object.")
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B))
    stop("B must be a positive integer.")

  states <- names(object)
  P <- object@transitionMatrix

  if (is.matrix(data)) {
    observed <- data
    if (is.null(rownames(observed)) || is.null(colnames(observed)))
      stop("A transition-count matrix must have row and column names.")
  } else if (is.atomic(data)) {
    observed <- createSequenceMatrix(stringchar = data, possibleStates = states)
  } else stop("data must be an empirical sequence or a transition-count matrix.")

  if (anyNA(observed) || any(observed < 0) || any(observed != floor(observed)))
    stop("Observed transition counts must be non-negative integers.")
  if (!setequal(rownames(observed), states) || !setequal(colnames(observed), states))
    stop("Empirical and theoretical transition matrices must have the same state support.")
  observed <- observed[states, states, drop = FALSE]

  row_totals <- rowSums(observed)
  expected <- sweep(P, 1L, row_totals, "*")
  informative <- row_totals > 0
  support <- outer(informative, rep(TRUE, length(states))) & (P > 0)

  if (any(observed > 0 & P == 0)) {
    result <- list(statistic = Inf, dof = NA_integer_, p.value = 0,
                   method = method, observed = observed,
                   expected = expected, transitionMatrix = P)
    if (verbose)
      cat("At least one observed transition has zero theoretical probability; the null hypothesis is rejected.\n")
    return(invisible(result))
  }

  statistic_fun <- function(obs) {
    positive <- (P > 0) & (obs > 0 | expected > 0)
    if (method == "Pearson")
      return(sum((obs[positive] - expected[positive])^2 / expected[positive]))
    2 * sum(obs[positive & obs > 0] *
              log(obs[positive & obs > 0] / expected[positive & obs > 0]))
  }
  statistic <- statistic_fun(observed)

  dof <- sum(vapply(which(informative),
                     function(i) sum(P[i, ] > 0) - 1L, integer(1)))

  if (method != "simulation") {
    p.value <- if (dof > 0) pchisq(statistic, df = dof, lower.tail = FALSE) else NA_real_
    result <- list(statistic = unname(statistic), dof = dof, p.value = p.value,
                   method = method, observed = observed, expected = expected,
                   transitionMatrix = P)
  } else {
    if (!is.null(seed)) set.seed(seed)
    simulated <- numeric(B)
    for (b in seq_len(B)) {
      sim <- matrix(0, nrow(P), ncol(P), dimnames = dimnames(P))
      for (i in which(informative))
        sim[i, ] <- as.vector(rmultinom(1L, row_totals[i], P[i, ]))
      expected_b <- expected
      # statistic_fun uses the observed expected table; construct the
      # conditional expected table once under the fixed theoretical P.
      positive_b <- (P > 0) & (sim > 0 | expected_b > 0)
      if (method == "Pearson") {
        simulated[b] <- sum((sim[positive_b] - expected_b[positive_b])^2 /
                            expected_b[positive_b])
      } else {
        pos <- positive_b & sim > 0
        simulated[b] <- 2 * sum(sim[pos] * log(sim[pos] / expected_b[pos]))
      }
    }
    p.value <- (1 + sum(simulated >= statistic)) / (B + 1)
    result <- list(statistic = unname(statistic), dof = dof, p.value = p.value,
                   method = method, observed = observed, expected = expected,
                   transitionMatrix = P, simulations = simulated, B = B)
  }

  if (verbose) cat("Testing compatibility of empirical and theoretical transition probabilities\n",
                   "Method:", method, "\nStatistic:", statistic,
                   "\nDegrees of freedom:", dof, "\np-value:", p.value, "\n")
  invisible(result)
}

#' Test homogeneity of transition probabilities across empirical sequences
#'
#' Tests whether several empirical transition matrices can be regarded as
#' realizations of the same time-homogeneous Markov chain.
#'
#' @rdname statisticalTests
#' @family statisticalTests
#' @param inputList A list whose elements are empirical sequences or matrices.
#' @param method Test statistic: `"G"`, `"Pearson"`, or `"simulation"`.
#' @param B Number of Monte Carlo replicates for `method = "simulation"`.
#' @param seed Optional random seed.
#' @param verbose Should results be printed?
#' @return A list containing statistic, degrees of freedom, p-value, method,
#'   pooled transition counts, and individual transition counts.
#' @references Kullback, S., Kupperman, M. and Ku, H. H. (1962). Tests for
#'   Contingency Tables and Markov Chains. *Technometrics*, 4(4), 573--608.
#' @export
verifyHomogeneity <- function(inputList,
                              method = c("G", "Pearson", "simulation"),
                              B = 9999, seed = NULL, verbose = TRUE) {
  method <- match.arg(method)
  if (!is.list(inputList) || length(inputList) < 2L)
    stop("inputList must contain at least two sequences or matrices.")
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B))
    stop("B must be a positive integer.")

  mats <- lapply(inputList, function(x) {
    if (is.matrix(x)) {
      if (is.null(rownames(x)) || is.null(colnames(x)))
        stop("Transition-count matrices must have row and column names.")
      m <- x
    } else if (is.atomic(x)) m <- createSequenceMatrix(stringchar = x)
    else stop("Each element must be a sequence or transition-count matrix.")
    if (anyNA(m) || any(m < 0) || any(m != floor(m)))
      stop("Transition counts must be non-negative integers.")
    m
  })

  states <- sort(unique(unlist(lapply(mats, rownames))))
  mats <- lapply(mats, function(m) {
    out <- matrix(0, length(states), length(states), dimnames = list(states, states))
    out[rownames(m), colnames(m)] <- m
    out
  })
  pooled <- Reduce("+", mats)
  statistic <- 0
  dof <- 0L

  for (i in seq_along(states)) {
    totals <- vapply(mats, function(m) sum(m[i, ]), numeric(1))
    active <- totals > 0
    if (sum(active) < 2L) next
    pooledRow <- pooled[i, ]
    destinations <- pooledRow > 0
    k <- sum(destinations)
    if (k < 2L) next
    p <- pooledRow[destinations] / sum(pooledRow)
    for (s in which(active)) {
      obs <- mats[[s]][i, destinations]
      exp <- totals[s] * p
      if (method == "G") {
        positive <- obs > 0
        statistic <- statistic + 2 * sum(obs[positive] * log(obs[positive] / exp[positive]))
      } else {
        statistic <- statistic + sum((obs - exp)^2 / exp)
      }
    }
    dof <- dof + (sum(active) - 1L) * (k - 1L)
  }

  if (method == "simulation") {
    if (!is.null(seed)) set.seed(seed)
    simulated <- numeric(B)
    for (b in seq_len(B)) {
      sim_mats <- lapply(mats, function(m) matrix(0, nrow(m), ncol(m),
                                                    dimnames = dimnames(m)))
      for (i in seq_along(states)) {
        pooledRow <- pooled[i, ]
        destinations <- pooledRow > 0
        total_pool <- sum(pooledRow)
        if (total_pool == 0) next
        p <- pooledRow[destinations] / total_pool
        for (s in seq_along(mats)) {
          total_s <- sum(mats[[s]][i, ])
          if (total_s > 0)
            sim_mats[[s]][i, destinations] <- as.vector(rmultinom(1L, total_s, p))
        }
      }

      stat_b <- 0
      for (i in seq_along(states)) {
        totals_b <- vapply(sim_mats, function(m) sum(m[i, ]), numeric(1))
        active_b <- totals_b > 0
        if (sum(active_b) < 2L) next
        pooled_b <- Reduce("+", sim_mats)[i, ]
        destinations_b <- pooled_b > 0
        k_b <- sum(destinations_b)
        if (k_b < 2L) next
        p_b <- pooled_b[destinations_b] / sum(pooled_b)
        for (s in which(active_b)) {
          obs_b <- sim_mats[[s]][i, destinations_b]
          exp_b <- totals_b[s] * p_b
          if (method == "Pearson") {
            stat_b <- stat_b + sum((obs_b - exp_b)^2 / exp_b)
          } else {
            pos <- obs_b > 0
            stat_b <- stat_b + 2 * sum(obs_b[pos] * log(obs_b[pos] / exp_b[pos]))
          }
        }
      }
      simulated[b] <- stat_b
    }
    p.value <- (1 + sum(simulated >= statistic)) / (B + 1)
  } else {
    p.value <- if (dof > 0) pchisq(statistic, df = dof, lower.tail = FALSE) else NA_real_
    simulated <- NULL
  }

  result <- list(statistic = unname(statistic), dof = dof, p.value = p.value,
                 method = method, pooled = pooled, transitionCounts = mats)
  if (method == "simulation") {
    result$simulations <- simulated
    result$B <- B
  }
  if (verbose) cat("Testing homogeneity of transition probabilities across sequences\n",
                   "Method:", method, "\nStatistic:", statistic,
                   "\nDegrees of freedom:", dof, "\np-value:", p.value, "\n")
  invisible(result)
}
