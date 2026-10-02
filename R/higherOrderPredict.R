#' Next-state probabilities and simulation for higher order Markov chains
#'
#' @description \code{higherOrderPredict} returns the distribution of the next
#'   state given the most recent states, and \code{higherOrderSimulate} draws a
#'   sequence of states, under a higher order model fitted by
#'   \code{\link{fitHigherOrder}} or \code{\link{fitMTD}}.
#'
#' @details Both functions use the transition probabilities of the fitted model
#'   of order \eqn{k},
#'   \deqn{P(x_t = j \mid x_{t-1}, \dots, x_{t-k}) = \sum_{i=1}^{k} \lambda_i\, Q_i[j, x_{t-i}],}
#'   the same that enter \code{\link{higherOrderLogLik}}. For a fit returned
#'   by \code{\link{fitMTD}} all the \eqn{Q_i} are the single MTD transition
#'   matrix. Only the last \eqn{k} states of a history are used, the last
#'   element being the most recent.
#'
#'   A model of order \eqn{k} needs \eqn{k} previous states, so \code{t0}
#'   (and every history) must contain at least \code{order} states. The
#'   probabilities are normalized to sum to one, which only matters when the
#'   weights of a least squares fit sum to one up to the optimizer's tolerance.
#'
#' @param fit The list returned by \code{\link{fitHigherOrder}} or
#'   \code{\link{fitMTD}}.
#' @param history The most recent states, oldest first: a vector of at least
#'   \code{order} states, or a matrix (or data frame) with one history per row.
#' @param n Number of states to simulate.
#' @param t0 The states preceding the simulated ones, oldest first; at least
#'   \code{order} states.
#' @param include.t0 Should \code{t0} be included at the beginning of the
#'   returned sequence?
#'
#' @return \code{higherOrderPredict}: a named vector of next-state
#'   probabilities for a single history, or a matrix with one row per history.
#'   \code{higherOrderSimulate}: a character vector of \code{n} states,
#'   preceded by \code{t0} if \code{include.t0 = TRUE}.
#'
#' @seealso \code{\link{fitHigherOrder}}, \code{\link{fitMTD}},
#'   \code{\link{higherOrderLogLik}}, \code{\link{rmarkovchain}}
#'
#' @examples
#' wind <- read.csv(system.file("extdata", "koeberg_wind.csv",
#'                              package = "markovchain"))$state
#' fit <- fitMTD(wind, order = 2)
#' # next wind direction after directions 1 and then 2
#' higherOrderPredict(fit, c(1, 2))
#' # several histories at once
#' higherOrderPredict(fit, rbind(c(1, 1), c(2, 2), c(4, 1)))
#' # simulate one day of hourly directions starting from the last two observed
#' set.seed(1)
#' higherOrderSimulate(24, fit, t0 = tail(wind, 2))
#'
#' # the same works for fitHigherOrder()
#' data(rain)
#' fit2 <- fitHigherOrder(rain$rain, order = 2, method = "mle")
#' higherOrderPredict(fit2, c("0", "6+"))
#' @name higherOrderPredict
#' @export
higherOrderPredict <- function(fit, history) {
  model <- .higherOrderModel(fit)
  if (is.data.frame(history)) history <- as.matrix(history)
  if (is.matrix(history)) {
    out <- t(apply(history, 1L, function(h) .higherOrderNext(model, h)))
    dimnames(out) <- list(rownames(history), model$states)
    return(out)
  }
  .higherOrderNext(model, history)
}

#' @rdname higherOrderPredict
#' @export
higherOrderSimulate <- function(n, fit, t0, include.t0 = FALSE) {
  if (length(n) != 1L || is.na(n) || n < 0 || n != floor(n)) {
    stop("n must be a non-negative integer")
  }
  model <- .higherOrderModel(fit)
  k <- length(model$lambda)
  history <- .higherOrderHistory(model, t0)
  out <- character(n)
  past <- utils::tail(history, k)
  for (t in seq_len(n)) {
    p <- .higherOrderProbabilities(model, past)
    out[t] <- model$states[sample.int(length(p), 1L, prob = p)]
    past <- c(past[-1L], out[t])
  }
  if (isTRUE(include.t0)) c(history, out) else out
}

# checks the fit and extracts what the predictions need
.higherOrderModel <- function(fit) {
  if (!is.list(fit) || is.null(fit$lambda) || !is.list(fit$Q) ||
      length(fit$Q) != length(fit$lambda) || length(fit$lambda) < 1L) {
    stop("fit must be the list returned by fitHigherOrder() or fitMTD()")
  }
  states <- rownames(fit$Q[[1L]])
  if (is.null(states)) stop("fit$Q must have state names")
  list(lambda = as.numeric(fit$lambda), Q = fit$Q, states = states)
}

.higherOrderHistory <- function(model, history) {
  if (!is.atomic(history) || anyNA(history)) {
    stop("a history must be a vector of states without missing values")
  }
  history <- as.character(history)
  k <- length(model$lambda)
  if (length(history) < k) {
    stop("a history must contain at least ", k, " states (the order of the model)")
  }
  if (!all(history %in% model$states)) stop("a history contains states unknown to the fitted model")
  history
}

# probabilities of the next state given the last k states (oldest first)
.higherOrderProbabilities <- function(model, past) {
  k <- length(model$lambda)
  p <- numeric(length(model$states))
  for (i in seq_len(k)) {
    p <- p + model$lambda[i] * model$Q[[i]][, past[k - i + 1L]]
  }
  p / sum(p)
}

.higherOrderNext <- function(model, history) {
  history <- .higherOrderHistory(model, history)
  p <- .higherOrderProbabilities(model, utils::tail(history, length(model$lambda)))
  setNames(as.numeric(p), model$states)
}
