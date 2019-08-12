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
#' @note This function is written in Rcpp.
#'
#' @examples
#' sequence<-c("a", "a", "b", "b", "a", "c", "b", "a", "b", "c", "a", "b",
#'             "c", "a", "b", "c", "a", "b", "a", "b")
#' fitHigherOrder(sequence)
#'
#' @export
fitHigherOrder<-function(sequence, order = 2) {
  # prbability of each states of sequence
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
  return(out)
}
