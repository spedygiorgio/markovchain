#' Test independence of consecutive states of an empirical sequence
#'
#' Tests the null hypothesis that consecutive observations are independent,
#' \eqn{P(X_{t+1} = j \mid X_t = i) = P(X_{t+1} = j)} for all \eqn{i, j},
#' against the alternative of a first-order Markov chain. This is the
#' classical test of Anderson and Goodman (1957): a chi-squared (or
#' likelihood-ratio) test of independence applied to the table of observed
#' one-step transition counts, whose rows are the state at time \eqn{t} and
#' whose columns are the state at time \eqn{t + 1}.
#'
#' Only states actually observed as a departure state (rows) or as an arrival
#' state (columns) contribute to the degrees of freedom, which are
#' \eqn{(r - 1)(c - 1)} for \eqn{r} such rows and \eqn{c} such columns. When
#' \eqn{r} or \eqn{c} is 1 (for instance a constant sequence) the test is not
#' defined: the degrees of freedom are 0 and the p-value is \code{NA}.
#'
#' The statistic is asymptotic and, as for any chi-squared test on a table of
#' counts, unreliable when many expected counts are small (say below 5).
#' Successive transitions overlap (each observation is the arrival state of one
#' transition and the departure state of the next); this is the standard
#' treatment of the Anderson-Goodman test and is asymptotically valid under
#' the null hypothesis.
#'
#' @family statisticalTests
#' @param sequence An empirical sequence of states (at least three
#'   observations, without missing values).
#' @param method Test statistic: \code{"Pearson"} (default) for the chi-squared
#'   statistic or \code{"G"} for the likelihood-ratio statistic.
#' @param verbose Should test results be printed?
#' @return An \code{htest} object, returned invisibly, with the additional
#'   components \code{observed} (transition counts, rows are departure states) and
#'   \code{expected} (counts expected under independence).
#' @references Anderson, T. W. and Goodman, L. A. (1957). Statistical inference
#' about Markov chains. \emph{The Annals of Mathematical Statistics}, 28(1), 89--110.
#' @seealso \code{\link{verifyMarkovProperty}}, \code{\link{assessOrder}},
#'   \code{\link{assessStationarity}}
#' @examples
#' # an independent sequence: the test should not reject
#' set.seed(1)
#' iid <- sample(c("a", "b", "c"), 500, replace = TRUE)
#' assessIndependence(iid)
#'
#' # a strongly dependent (Markov) sequence: the test rejects
#' mc <- new("markovchain", states = c("a", "b"),
#'           transitionMatrix = matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE))
#' dep <- rmarkovchain(500, mc, t0 = "a")
#' assessIndependence(dep, method = "G")
#' @export
assessIndependence <- function(sequence, method = c("Pearson", "G"), verbose = TRUE) {
  method <- match.arg(method)
  if (length(sequence) < 3L) stop("sequence must contain at least three observations.")
  if (anyNA(sequence)) stop("sequence must not contain missing values.")
  data.name <- deparse(substitute(sequence))
  counts <- createSequenceMatrix(as.character(sequence))
  row_totals <- rowSums(counts)
  col_totals <- colSums(counts)
  observed <- counts[row_totals > 0, col_totals > 0, drop = FALSE]
  n <- sum(observed)
  expected <- outer(row_totals[row_totals > 0], col_totals[col_totals > 0]) / n
  dimnames(expected) <- dimnames(observed)
  dof <- (nrow(observed) - 1L) * (ncol(observed) - 1L)
  if (method == "Pearson") {
    statistic <- sum((observed - expected)^2 / expected)
  } else {
    positive <- observed > 0
    statistic <- 2 * sum(observed[positive] * log(observed[positive] / expected[positive]))
  }
  result <- list(statistic = statistic, dof = dof,
                 p.value = if (dof > 0L) pchisq(statistic, dof, lower.tail = FALSE) else NA_real_,
                 observed = observed, expected = expected)
  result <- .asHtest(result, method, data.name)
  result$method <- paste(result$method, "for independence of consecutive states")
  if (verbose) print(result)
  invisible(result)
}
