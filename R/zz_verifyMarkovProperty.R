#' Test the first-order Markov property of an empirical sequence
#'
#' Tests the null hypothesis that the conditional distribution of the next
#' state depends only on the current state, against a second-order Markov
#' alternative. Equivalently, the test assesses
#' `X[t-1] independent of X[t+1] | X[t]`.
#'
#' @rdname statisticalTests
#' @family statisticalTests
#' @param sequence An empirical sequence of states.
#' @param method Test statistic: `"G"`, `"Pearson"`, or `"simulation"`.
#' @param B Number of Monte Carlo replicates for `method = "simulation"`.
#' @param seed Optional random seed.
#' @param verbose Should results be printed?
#'
#' @return A list with statistic, degrees of freedom, p-value, observed and
#' expected triple-transition counts, and method. For simulation, simulated
#' statistics are also returned.
#'
#' @references
#' Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about
#' Markov chains. Annals of Mathematical Statistics, 28(1), 89--110.
#'
#' @export
verifyMarkovProperty <- function(sequence,
                                 method = c("G", "Pearson", "simulation"),
                                 B = 9999,
                                 seed = NULL,
                                 verbose = TRUE) {
  method <- match.arg(method)

  if (length(sequence) < 4L)
    stop("sequence must contain at least four observations.")
  if (anyNA(sequence))
    stop("sequence must not contain missing values.")
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B))
    stop("B must be a positive integer.")

  states <- unique(sequence)
  r <- length(states)
  n <- length(sequence)
  idx <- match(sequence, states)

  Nijk <- array(0L, dim = c(r, r, r),
                dimnames = list(states, states, states))
  linear <- idx[-c(1L, 2L)] +
    (idx[2:(n - 1L)] - 1L) * r +
    (idx[1:(n - 2L)] - 1L) * r * r
  Nijk[] <- tabulate(linear, nbins = r^3L)

  Nij <- apply(Nijk, c(1L, 2L), sum)
  Njk <- apply(Nijk, c(2L, 3L), sum)
  Nj <- rowSums(Njk)

  expected <- array(0, dim = c(r, r, r),
                    dimnames = dimnames(Nijk))
  for (j in seq_len(r)) {
    if (Nj[j] > 0L)
      expected[, j, ] <- outer(Nij[, j], Njk[j, ]) / Nj[j]
  }

  statistic_fun <- function(observed, expected, method) {
    positive_expected <- expected > 0
    if (method == "Pearson") {
      return(sum((observed[positive_expected] -
                  expected[positive_expected])^2 /
                 expected[positive_expected]))
    }
    positive_observed <- observed > 0
    return(2 * sum(observed[positive_observed] *
                   log(observed[positive_observed] /
                       expected[positive_observed])))
  }

  statistic <- statistic_fun(
    Nijk, expected,
    method = if (method == "Pearson") "Pearson" else "G"
  )

  # Conditional-independence degrees of freedom. For each current state j,
  # the contribution is (number of observed previous states - 1) times
  # (number of observed future states - 1). For complete support this is
  # r * (r - 1)^2.
  df <- 0L
  for (j in seq_len(r)) {
    previous <- sum(Nij[, j] > 0)
    future <- sum(Njk[j, ] > 0)
    if (previous > 0L && future > 0L)
      df <- df + (previous - 1L) * (future - 1L)
  }

  transition_matrix <- matrix(
    0, r, r, dimnames = list(states, states)
  )
  informative <- Nj > 0
  transition_matrix[informative, ] <- sweep(
    Njk[informative, , drop = FALSE], 1L, Nj[informative], "/"
  )

  if (method != "simulation") {
    p.value <- pchisq(statistic, df = df, lower.tail = FALSE)
    result <- list(
      statistic = statistic,
      dof = df,
      p.value = p.value,
      method = method,
      observed = Nijk,
      expected = expected,
      transitionMatrix = transition_matrix
    )
  } else {
    if (!all(informative)) {
      stop(
        "method = 'simulation' requires at least one observed outgoing " ,
        "transition from every state appearing in the sequence."
      )
    }

    if (!is.null(seed)) set.seed(seed)

    # Parametric bootstrap under the fitted first-order Markov model.
    # The transition matrix is re-estimated in every replicate, so parameter
    # estimation is included in the simulated null distribution.
    simulated <- numeric(B)
    initial <- idx[1L]
    P <- transition_matrix

    for (b in seq_len(B)) {
      x <- integer(n)
      x[1L] <- initial
      for (t in seq_len(n - 1L))
        x[t + 1L] <- sample.int(r, 1L, prob = P[x[t], ])

      a <- x[-c(1L, 2L)] +
        (x[2:(n - 1L)] - 1L) * r +
        (x[1:(n - 2L)] - 1L) * r * r
      Ns <- array(tabulate(a, nbins = r^3L), dim = c(r, r, r))
      Nijs <- apply(Ns, c(1L, 2L), sum)
      Njks <- apply(Ns, c(2L, 3L), sum)
      Njs <- rowSums(Njks)
      Es <- array(0, dim = c(r, r, r))
      for (j in seq_len(r)) {
        if (Njs[j] > 0L)
          Es[, j, ] <- outer(Nijs[, j], Njks[j, ]) / Njs[j]
      }

      simulated[b] <- statistic_fun(
        Ns, Es,
        method = if (method == "Pearson") "Pearson" else "G"
      )
    }

    p.value <- (1 + sum(simulated >= statistic)) / (B + 1)
    result <- list(
      statistic = statistic,
      dof = df,
      p.value = p.value,
      method = method,
      observed = Nijk,
      expected = expected,
      transitionMatrix = transition_matrix,
      simulations = simulated,
      B = B,
      initial.state = states[initial]
    )
  }

  if (verbose) {
    cat("Testing compatibility with a first-order Markov model\n")
    cat("Method:", method, "\n")
    cat("Statistic:", statistic, "\n")
    cat("Degrees of freedom:", df, "\n")
    cat("p-value:", p.value, "\n")
  }

  invisible(result)
}