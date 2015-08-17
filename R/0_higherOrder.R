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

fn1=function(x)
{
  exp(x[1]*x[2]*x[3]*x[4]*x[5])
}
eqn1=function(x){
  z1=x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[4]*x[4]+x[5]*x[5]
  z2=x[2]*x[3]-5*x[4]*x[5]
  z3=x[1]*x[1]*x[1]+x[2]*x[2]*x[2]
  return(c(z1,z2,z3))
}

.seq2freqProb<-function(sequence) {
  n<-length(sequence)
  states<-unique(sequence)
  nstates<-length(states)
  v<-vector(mode="numeric", length=nstates)
  names(v)<-states
  # print(v)
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

# fitHigherOrder<-function(frequency, sequencelist, order) {
fitHigherOrder<-function(sequence, order = 2) {
  # print(sequence)
  F<-list()
  Q<-list()
  for(o in 1:order) {
    F[[o]]<-.seq2matHigh(sequence, o)
    # print(colSums(F[[o]]))
    Q[[o]]<-sweep(F[[o]], 2, colSums(F[[o]]), '/')
  }
  # print(F)
  # print(Q)
  X<-.seq2freqProb(sequence)
  # print(X)
  x0 <- c(-2, 2, 2, -1, -1)
  model<-Rsolnp::solnp(x0, fun = fn1, eqfun = eqn1, eqB = c(10, 0, -1), control=list(trace=0))
  # print(model)
  
  # X=as.numeric(frequency/sum(frequency))
  # ll=vector()
  # Q=alply(.data=seq(1,order,1), .margins=1, .fun=.getQ, sequencelist)
  # transitions=llply(.data=Q, .fun=function(q) q$transition)
  # ll=laply(.data=Q, .fun=function(q) q$ll)
  # QX=llply(.data=transitions, .fun=function(tr) as.matrix(tr)%*%X) 
#   environment(.foo)=environment()
#   params=rep(1/order, order)
#   model=Rsolnp::solnp(pars=params, fun=.foo, eqfun=.constr, eqB=1, LB=rep(0, order), control=list(trace=0))
  lambda=model$pars
  return(lambda)
}
