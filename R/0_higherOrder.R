# define higher order Markov Chain class

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

# fit higher order markov chain
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
