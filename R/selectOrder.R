#' Select the order of a Markov chain by information criteria
#'
#' @description Fits by maximum likelihood fully parameterized Markov chains of
#'   order \eqn{0, 1, \dots,} \code{maxOrder} to an empirical sequence (or to a
#'   list of independent sequences), all on the same observations, and selects
#'   the order that minimizes the BIC or the AIC (Tong, 1975; Katz, 1981).
#'
#' @details A Markov chain of order \eqn{k} on \eqn{r} states has a transition
#'   probability for each context (the \eqn{k} previous states) and next state;
#'   order 0 is independence. Its maximum likelihood estimates are the observed
#'   transition frequencies, and its log-likelihood is
#'   \eqn{\sum N(c, j) \log\{N(c, j) / N(c)\}}, where \eqn{N(c, j)} counts
#'   context \eqn{c} followed by state \eqn{j}.
#'
#'   Information criteria are only comparable when every model is evaluated on
#'   the same observations. An order-\eqn{k} chain can only predict an
#'   observation from position \eqn{k + 1} onwards, so all orders are
#'   evaluated on the observations from \code{start} (by default
#'   \code{maxOrder + 1}) to the end of each sequence, the earlier ones being
#'   only used as contexts. Berchtold and Raftery (2002), for instance, condition
#'   on the first 14 observations (\code{start = 15}). For a list of sequences
#'   the counts are pooled, the observations before \code{start} of each
#'   sequence are only used as contexts, and no context crosses from one
#'   sequence to the next.
#'
#'   With \code{parameters = "full"} (the default) an order-\eqn{k} chain has
#'   \eqn{r^k (r - 1)} free parameters, as in Tong (1975) and Katz (1981).
#'   With \code{parameters = "observed"} only the probabilities that are not
#'   estimated as zero are counted, that is the number of distinct next states
#'   minus one summed over the observed contexts; this is the convention of
#'   Berchtold and Raftery (2002), which is less penalizing when many
#'   transitions are never observed. In both cases BIC uses the number of
#'   observations entering the likelihood.
#'
#'   The table also reports the likelihood-ratio statistic of each order
#'   against the previous one, \eqn{G = 2 (\ell_k - \ell_{k-1})}, with degrees
#'   of freedom equal to the difference in the number of parameters and an
#'   asymptotic chi-squared p-value (Anderson and Goodman, 1957). These tests
#'   are not adjusted for multiplicity, and they and the criteria become
#'   unreliable when the number of contexts \eqn{r^k} is not small compared with
#'   the number of observations: BIC is consistent for the order (Csiszar and
#'   Shields, 2000), whereas AIC tends to select too high an order in long
#'   sequences (Katz, 1981).
#'
#' @param sequence An empirical sequence of states (a vector coercible to
#'   character, without missing values), or a list of such sequences.
#' @param maxOrder The highest order considered, a non-negative integer.
#' @param criterion The criterion used to select the order: \code{"BIC"}
#'   (default) or \code{"AIC"}.
#' @param start Index of the first observation of each sequence entering the
#'   likelihood, at least \code{maxOrder + 1} (the default).
#' @param parameters How the free parameters are counted: \code{"full"}
#'   (default) or \code{"observed"}, see Details.
#'
#' @return A list with components
#'   \item{order}{the order selected by \code{criterion}}
#'   \item{criterion}{the criterion used}
#'   \item{table}{a data frame with one row per order: \code{order},
#'     \code{logLik}, \code{npar}, \code{AIC}, \code{BIC}, and the
#'     likelihood-ratio test against the previous order (\code{LR},
#'     \code{df}, \code{p.value}; \code{NA} for order 0)}
#'   \item{nobs}{number of observations entering the likelihood}
#'   \item{start, parameters}{as used}
#'
#' @references
#' Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about Markov
#' chains. The Annals of Mathematical Statistics, 28(1), 89-110.
#'
#' Tong, H. (1975). Determination of the order of a Markov chain by Akaike's
#' information criterion. Journal of Applied Probability, 12(3), 488-497.
#'
#' Katz, R. W. (1981). On some criteria for estimating the order of a Markov
#' chain. Technometrics, 23(3), 243-249.
#'
#' Csiszar, I. and Shields, P. C. (2000). The consistency of the BIC Markov
#' order estimator. The Annals of Statistics, 28(6), 1601-1619.
#'
#' Berchtold, A. and Raftery, A. E. (2002). The mixture transition distribution
#' model for high-order Markov chains and non-Gaussian time series.
#' Statistical Science, 17(3), 328-356.
#'
#' @seealso \code{\link{assessOrder}}, \code{\link{verifyMarkovProperty}},
#'   \code{\link{fitHigherOrder}}, \code{\link{fitMTD}},
#'   \code{\link{higherOrderLogLik}}
#'
#' @examples
#' # Alofi rainfall: three states
#' data(rain)
#' selectOrder(rain$rain, maxOrder = 3)$table
#'
#' # Koeberg wind directions with the conventions of Berchtold and Raftery (2002)
#' wind <- read.csv(system.file("extdata", "koeberg_wind.csv",
#'                              package = "markovchain"))$state
#' sel <- selectOrder(wind, maxOrder = 3, start = 15, parameters = "observed")
#' sel$order
#' sel$table
#' @export
selectOrder <- function(sequence, maxOrder = 3, criterion = c("BIC", "AIC"),
                        start = NULL, parameters = c("full", "observed")) {
  criterion <- match.arg(criterion)
  parameters <- match.arg(parameters)
  sequences <- if (is.list(sequence)) sequence else list(sequence)
  if (length(sequences) == 0L) stop("sequence must contain at least one sequence")
  for (s in sequences) {
    if (!is.atomic(s) || anyNA(s)) {
      stop("each sequence must be a vector of states without missing values")
    }
  }
  sequences <- lapply(sequences, as.character)
  if (length(maxOrder) != 1L || is.na(maxOrder) || !is.finite(maxOrder) ||
      maxOrder < 0 || maxOrder != floor(maxOrder)) {
    stop("maxOrder must be a non-negative integer")
  }
  maxOrder <- as.integer(maxOrder)
  if (is.null(start)) start <- maxOrder + 1L
  if (length(start) != 1L || is.na(start) || !is.finite(start) ||
      start != floor(start) || start < maxOrder + 1L) {
    stop("start must be an integer of at least maxOrder + 1")
  }
  start <- as.integer(start)
  lengths <- vapply(sequences, length, integer(1))
  sequences <- sequences[lengths >= start]
  if (length(sequences) == 0L) {
    stop("no sequence has observations from position start onwards")
  }
  r <- length(unique(unlist(sequences)))
  if (r < 2L) stop("the sequences must contain at least two different states")

  # next states and, for each order, the context of every observation used
  nextState <- unlist(lapply(sequences, function(s) s[start:length(s)]))
  nobs <- length(nextState)
  context <- function(k) {
    if (k == 0L) return(rep("", nobs))
    unlist(lapply(sequences, function(s) {
      times <- start:length(s)
      lagged <- lapply(k:1, function(g) s[times - g])
      do.call(paste, c(lagged, sep = "\r"))
    }))
  }

  orders <- 0:maxOrder
  logLik <- npar <- numeric(length(orders))
  for (k in orders) {
    ctx <- context(k)
    counts <- table(ctx, nextState)
    contextTotals <- rowSums(counts)
    positive <- counts > 0
    logLik[k + 1L] <- sum(counts[positive] *
                            log(counts[positive] / (contextTotals %o% rep(1, ncol(counts)))[positive]))
    npar[k + 1L] <- if (parameters == "full") r^k * (r - 1) else sum(rowSums(positive) - 1)
  }
  if (r^maxOrder > nobs) {
    warning("the number of contexts of the highest order (", r^maxOrder,
            ") exceeds the number of observations (", nobs,
            "): the criteria and tests for the highest orders are unreliable")
  }
  df <- c(NA, diff(npar))
  LR <- c(NA, 2 * diff(logLik))
  p.value <- ifelse(!is.na(df) & df > 0, stats::pchisq(pmax(LR, 0), df, lower.tail = FALSE), NA_real_)
  table <- data.frame(order = orders, logLik = logLik, npar = npar,
                      AIC = -2 * logLik + 2 * npar,
                      BIC = -2 * logLik + log(nobs) * npar,
                      LR = LR, df = df, p.value = p.value)
  list(order = orders[which.min(table[[criterion]])], criterion = criterion,
       table = table, nobs = nobs, start = start, parameters = parameters)
}
