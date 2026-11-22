#' Fit a continuous time Markov chain
#'
#' Maximum likelihood estimation of a continuous time Markov
#' chain from one observed trajectory: the sequence of states visited and
#' the times at which they were entered.
#'
#' @param data A list of two vectors of equal length: the states visited
#'   (character) and the times at which they were entered (numeric,
#'   increasing, the first equal to 0).
#' @param byrow Logical: return the generator, and the transition matrix of
#'   the embedded jump chain in the confidence intervals, with the rates out
#'   of each state in its row (default) or in its column?
#' @param name Optional name of the fitted chain.
#' @param confidencelevel Confidence level of the intervals (default 0.95).
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{estimate}}{the fitted \code{\link[=ctmc-class]{ctmc}} object;}
#'     \item{\code{errors}}{a list with \code{dtmcConfidenceInterval}, the confidence
#'       level and the lower and upper endpoint matrices of the transition
#'       probabilities of the embedded jump chain (as in \code{\link[=markovchainFit]{markovchainFit()}}),
#'       and \code{lambdaConfidenceInterval}, the lower and upper endpoints of the
#'       exit rate of each state.}
#'   }
#' @details The exit rate of state \eqn{i} is estimated as the number of
#'   departures from \eqn{i} divided by the total time spent in \eqn{i}, and
#'   the jump probabilities by the transition frequencies of the sequence of
#'   states; the generator is their product. The interval for each rate is the
#'   normal approximation \eqn{\hat\lambda_i (1 \pm z / \sqrt{n_i})},
#'   with \eqn{n_i} the number of departures from \eqn{i}, truncated at 0.
#' @references Yalamanchi, S. B. and Spedicato, G. A. (2015). Continuous time
#'   Markov chains. Vignette of the markovchain package.
#' @author Sai Bhargav Yalamanchi
#' @family continuous time Markov chains
#' @examples
#' data <- list(c("a", "b", "c", "a", "b", "a", "c", "b", "c"),
#'              c(0, 0.8, 2.1, 2.4, 4, 5, 5.9, 8.2, 9))
#' fit <- ctmcFit(data)
#' fit$estimate
#' fit$errors$lambdaConfidenceInterval
#' @export
ctmcFit <- function(data, byrow = TRUE, name = "", confidencelevel = 0.95) {
  .Call(`_markovchain_ctmcFit`, data, byrow, name, confidencelevel)
}
