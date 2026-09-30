#' @title Higher order Markov Chains class
#' @name HigherOrderMarkovChain-class
#' @description The S4 class that describes \code{HigherOrderMarkovChain} objects.
#' 
#' @export
setClass("HigherOrderMarkovChain", #class name
         representation(
           states = "character", 
           order = "numeric",
           transitions = "list", 
           name = "character"
         )
#          , prototype(states = c("a","b"), byrow = TRUE, # prototypizing
#                    transitionMatrix=matrix(data = c(0,1,1,0),
#                                            nrow=2, byrow=TRUE, dimnames=list(c("a","b"), c("a","b"))),
#                    name="Unnamed Markov chain")
)

# objective function to pass to solnp
.fn1=function(params)
{
  QX <- get("QX")
  X <- get("X")    
  error <- 0
  for (i in 1:length(QX)) {
    error <- error+(params[i] * QX[[i]]-X)
  }
  return(sum(error^2))
}

# equality constraint function to pass to solnp
.eqn1=function(params){
  return(sum(params))
}

#' @name fitHigherOrder
#' @aliases seq2freqProb seq2matHigh 
#' @title Functions to fit a higher order Markov chain
#'
#' @description Given a sequence of states arising from a stationary state, it
#'   fits the underlying Markov chain distribution with higher order.
#' @usage  
#' fitHigherOrder(sequence, order = 2)
#' seq2freqProb(sequence)
#' seq2matHigh(sequence, order)
#'
#' @param sequence A character list.
#' @param order Markov chain order
#' @return A list containing lambda, Q, and X.
#'
#' @references 
#' Ching, W. K., Huang, X., Ng, M. K., & Siu, T. K. (2013). Higher-order markov 
#' chains. In Markov Chains (pp. 141-176). Springer US.
#' 
#' Ching, W. K., Ng, M. K., & Fung, E. S. (2008). Higher-order multivariate
#' Markov chains and their applications. Linear Algebra and its Applications,
#' 428(2), 492-507.
#'
#' @author Giorgio Spedicato, Tae Seung Kang

#'
#' @examples
#' sequence<-c("a", "a", "b", "b", "a", "c", "b", "a", "b", "c", "a", "b",
#'             "c", "a", "b", "c", "a", "b", "a", "b")
#' fitHigherOrder(sequence)
#'
#' @export
fitHigherOrder<-function(sequence, order = 2) {
  if (!is.character(sequence) || length(sequence) < 2L || anyNA(sequence)) {
    stop("sequence must be a non-empty character vector without missing values")
  }
  if (length(order) != 1L || is.na(order) || !is.finite(order) ||
      order < 1 || order != floor(order) || order >= length(sequence)) {
    stop("order must be a positive integer smaller than the sequence length")
  }
  order <- as.integer(order)
  # prbability of each states of sequence
  if (requireNamespace("Rsolnp", quietly = TRUE)) {
  X <- seq2freqProb(sequence)
  
  # store h step transition matrix
  Q <- list()
  QX <- list()
  for(o in 1:order) {
    Q[[o]] <- seq2matHigh(sequence, o)
    QX[[o]] <- Q[[o]]%*%X
  }
  environment(.fn1) <- environment()
  params <- rep(1/order, order)
  model <- Rsolnp::solnp(params, fun=.fn1, eqfun=.eqn1, eqB=1, 
                         LB=rep(0, order), control=list(trace=0))
  lambda <- model$pars
  out <- list(lambda=lambda, Q=Q, X=X)
  } else {
    print("package Rsolnp unavailable")
    out <- NULL
  }
  return(out)
}


#' Log-likelihood, deviance and information criteria of a higher order Markov chain
#'
#' @description Evaluates the log-likelihood of an empirical sequence under the
#'   higher order Markov chain returned by \code{\link{fitHigherOrder}}, and
#'   derives the deviance, AIC and BIC, so that models of different orders can
#'   be compared.
#'
#' @details The fitted model is the mixture-transition-distribution model of
#'   Raftery (1985) in the form used by Ching et al.: the probability of moving
#'   to state \eqn{x_t} given the past is
#'   \deqn{P(x_t \mid x_{t-1}, \dots, x_{t-k}) = \sum_{i=1}^{k} \lambda_i\, Q_i[x_t, x_{t-i}],}
#'   where \eqn{Q_i} is the lag-\eqn{i} transition matrix (see
#'   \code{\link{seq2matHigh}}) and \eqn{k} is the order. The log-likelihood is
#'   the sum of the logarithms of these probabilities over the observations
#'   \eqn{t = } \code{start}, \eqn{\dots}, \eqn{T}, and the deviance is
#'   \eqn{-2} times the log-likelihood.
#'
#'   Two points matter when interpreting the output. First,
#'   \code{fitHigherOrder} chooses \eqn{\lambda} by least squares on the
#'   stationary distribution, not by maximum likelihood, so the value returned
#'   is the log-likelihood \emph{of the fitted model}, not the maximum
#'   attainable one. Second, a model of order \eqn{k} can only be evaluated from
#'   observation \eqn{k + 1} onwards; to compare orders on exactly the same data
#'   set \code{start} to \eqn{1 +} the largest order compared, otherwise the
#'   models are evaluated on different numbers of observations and neither the
#'   log-likelihood nor the information criteria are comparable.
#'
#'   The number of parameters used for AIC and BIC is
#'   \eqn{k\, r (r - 1) + (k - 1)}, that is \eqn{r (r - 1)} free probabilities
#'   for each of the \eqn{k} lag matrices plus the \eqn{k - 1} free weights,
#'   with \eqn{r} the number of states.
#'
#' @param sequence A character vector, the empirical sequence of states.
#' @param fit The list returned by \code{\link{fitHigherOrder}} for
#'   \code{sequence}. If \code{NULL}, \code{fitHigherOrder(sequence, order)} is
#'   computed.
#' @param order Order of the model to fit when \code{fit} is \code{NULL}
#'   (ignored otherwise; the order is then \code{length(fit$lambda)}).
#' @param start Index of the first observation included in the likelihood.
#'   Defaults to \code{order + 1}, the earliest observation a model of that
#'   order can predict.
#'
#' @return A list with components \code{logLik}, \code{deviance}, \code{AIC},
#'   \code{BIC}, \code{nobs} (number of observations entering the likelihood),
#'   \code{npar}, \code{order} and \code{start}. The log-likelihood is
#'   \code{-Inf} if an observed transition has probability zero under the model
#'   (possible only for a sequence other than the one the model was fitted on).
#'
#' @references
#' Raftery, A. E. (1985). A model for high-order Markov chains. Journal of the
#' Royal Statistical Society, Series B, 47(3), 528-539.
#'
#' Ching, W. K., Huang, X., Ng, M. K., & Siu, T. K. (2013). Higher-order markov
#' chains. In Markov Chains (pp. 141-176). Springer US.
#'
#' @seealso \code{\link{fitHigherOrder}}
#'
#' @examples
#' sequence <- c("a", "a", "b", "b", "a", "c", "b", "a", "b", "c", "a", "b",
#'               "c", "a", "b", "c", "a", "b", "a", "b")
#' # compare orders 1 and 2 on the same observations (start = 3)
#' if (requireNamespace("Rsolnp", quietly = TRUE)) {
#'   fit1 <- fitHigherOrder(sequence, order = 1)
#'   fit2 <- fitHigherOrder(sequence, order = 2)
#'   sapply(list(order1 = fit1, order2 = fit2), function(f)
#'     unlist(higherOrderLogLik(sequence, f, start = 3)[c("logLik", "deviance", "AIC", "BIC")]))
#' }
#'
#' @export
higherOrderLogLik <- function(sequence, fit = NULL, order = 2, start = NULL) {
  if (!is.character(sequence) || length(sequence) < 2L || anyNA(sequence)) {
    stop("sequence must be a non-empty character vector without missing values")
  }
  if (is.null(fit)) {
    fit <- fitHigherOrder(sequence, order)
    if (is.null(fit)) stop("package Rsolnp is required to fit the model")
  }
  if (!is.list(fit) || is.null(fit$lambda) || !is.list(fit$Q) ||
      length(fit$Q) != length(fit$lambda)) {
    stop("fit must be the list returned by fitHigherOrder()")
  }
  lambda <- as.numeric(fit$lambda)
  k <- length(lambda)
  states <- rownames(fit$Q[[1L]])
  if (is.null(states)) stop("fit$Q must have state names")
  idx <- match(sequence, states)
  if (anyNA(idx)) stop("sequence contains states unknown to the fitted model")
  n <- length(sequence)
  if (is.null(start)) start <- k + 1L
  if (length(start) != 1L || !is.finite(start) || start != floor(start) ||
      start < k + 1L || start > n) {
    stop("start must be an integer between order + 1 and the sequence length")
  }
  times <- start:n
  p <- numeric(length(times))
  for (i in seq_len(k)) {
    p <- p + lambda[i] * fit$Q[[i]][cbind(idx[times], idx[times - i])]
  }
  logLik <- if (any(p <= 0)) -Inf else sum(log(p))
  r <- length(states)
  npar <- k * r * (r - 1) + (k - 1)
  list(logLik = logLik, deviance = -2 * logLik,
       AIC = -2 * logLik + 2 * npar, BIC = -2 * logLik + log(length(times)) * npar,
       nobs = length(times), npar = npar, order = k, start = as.integer(start))
}
