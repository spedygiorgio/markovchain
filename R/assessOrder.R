#' Assess the Markov order of an empirical sequence
#'
#' Tests whether the conditional distribution of the future state depends on
#' the previous state after conditioning on the current state.
#'
#' @rdname statisticalTests
#' @family statisticalTests
#' @param sequence An empirical sequence of states.
#' @param verbose Should results be printed?
#' @export
assessOrder <- function(sequence, verbose = TRUE) {
  warning("The accuracy of the statistical inference functions has been questioned. It will be thoroughly investigated in future versions of the package.")
  n <- length(sequence)
  states <- unique(sequence)
  nelements <- length(states)
  TStat <- 0
  for (present in states) {
    mat <- matrix(0, nelements, nelements,
                  dimnames = list(states, states))
    for (i in 1:(n - 2)) {
      if (present == sequence[i + 1]) {
        past <- sequence[i]
        future <- sequence[i + 2]
        mat[past, future] <- mat[past, future] + 1
      }
    }
    res <- stats::chisq.test(mat)
    TStat <- TStat + res$statistic
  }
  k <- nelements
  df <- k * (k - 1)^2
  pvalue <- 1 - stats::pchisq(q = TStat, df)
  out <- list(statistic = TStat[[1]], p.value = pvalue[[1]])
  if (verbose) {
    cat("The assessOrder test statistic is: ", TStat, "\n")
    cat("The Chi-Square d.f. are: ", df, "\n")
    cat("The p-value is: ", pvalue, "\n")
  }
  invisible(out)
}
