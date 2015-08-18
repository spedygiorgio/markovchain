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

.seq2freqProb<-function(sequence) {
  n<-length(sequence)
  states<-unique(sequence)
  nstates<-length(states)
  v<-vector(mode="numeric", length=nstates)
  names(v)<-states
  for(i in 1:n) {
    v[sequence[i]]<-v[sequence[i]] + 1
  }
  return (v/sum(v))
}

.seq2matHigh<-function(sequence, order) {
  n<-length(sequence)
  states<-unique(sequence)
  nstates<-length(states)
  mat<-zeros(nstates)
  dimnames(mat)<-list(states, states)
  for(i in 1:(n-order)) {
    from<-sequence[i]
    to<-sequence[i+order]
    mat[to,from]<-mat[to,from]+1
  }
  return (mat)
}

.fn1=function(params)
{
  QX=get("QX")
  X=get("X")    
  error=0
  for (i in 1:length(QX)) {
    error=error+(params[i]*QX[[i]]-X)
  }
  return(sum(error^2))
}

.eqn1=function(params){
  return(sum(params))
}

fitHigherOrder<-function(sequence, order = 2) {
  X<-.seq2freqProb(sequence)
  F<-list()
  Q<-list()
  QX<-list()
  for(o in 1:order) {
    F[[o]]<-.seq2matHigh(sequence, o)
    Q[[o]]<-sweep(F[[o]], 2, colSums(F[[o]]), '/')
    QX[[o]]<-Q[[o]]%*%X
    # print(QX[[o]])
  }
  environment(.fn1)=environment()
  params<-rep(1/order, order)
  model<-Rsolnp::solnp(params, fun = .fn1, eqfun = .eqn1, eqB=1, LB=rep(0, order), control=list(trace=0))
  # model<-Rsolnp::solnp(x0, fun = .fn1, eqfun = .eqn1, eqB = c(10, 0, -1), control=list(trace=0))
  # print(model)
#   model=Rsolnp::solnp(pars=params, fun=.foo, eqfun=.constr, eqB=1, LB=rep(0, order), control=list(trace=0))
  lambda=model$pars
  return(lambda)
}
