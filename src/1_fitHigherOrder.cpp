#include <Rcpp.h>
using namespace Rcpp;

NumericVector seq2freqProb (CharacterVector sequence) {
  int n = sequence.size(); 
  CharacterVector states = unique(sequence).sort();
  int nstates = states.length();
  NumericVector v(nstates);
  v.names() = states;
  for(int i = 0; i < n; i ++) {
    v[std::string(sequence[i])] = v[std::string(sequence[i])] + 1;
  }
  NumericVector out = v/sum(v);
  out.names() = v.names();
  return out;
}

// 
// .seq2matHigh<-function(sequence, order) {
// n<-length(sequence)
//   states<-unique(sequence)
//   nstates<-length(states)
//   mat<-zeros(nstates)
//   dimnames(mat)<-list(states, states)
//   for(i in 1:(n-order)) {
//     from<-sequence[i]
//     to<-sequence[i+order]
//     mat[to,from]<-mat[to,from]+1
//   }
//   return (mat)
// }
// 
// .fn1=function(params)
// {
//   QX=get("QX")
//   X=get("X")    
//   error=0
//   for (i in 1:length(QX)) {
//     error=error+(params[i]*QX[[i]]-X)
//   }
//   return(sum(error^2))
// }
// 
// .eqn1=function(params){
// return(sum(params))
// }
// 
// fitHigherOrder<-function(sequence, order = 2) {
//   X<-.seq2freqProb(sequence)
//   F<-list()
//   Q<-list()
//   QX<-list()
//   for(o in 1:order) {
//     F[[o]]<-.seq2matHigh(sequence, o)
//     Q[[o]]<-sweep(F[[o]], 2, colSums(F[[o]]), '/')
//     QX[[o]]<-Q[[o]]%*%X
// # print(QX[[o]])
//   }
//   environment(.fn1)=environment()
//     params<-rep(1/order, order)
//     model<-Rsolnp::solnp(params, fun = .fn1, eqfun = .eqn1, eqB=1, LB=rep(0, order), control=list(trace=0))
// # model<-Rsolnp::solnp(x0, fun = .fn1, eqfun = .eqn1, eqB = c(10, 0, -1), control=list(trace=0))
// # print(model)
// #   model=Rsolnp::solnp(pars=params, fun=.foo, eqfun=.constr, eqB=1, LB=rep(0, order), control=list(trace=0))
//     lambda=model$pars
//     return(lambda)
// }

// [[Rcpp::export]]
void fitHigherOrderRcpp(SEXP sequence, int order = 2) {
  Rcout << "fitHigherOrder " << order << std::endl;
  // Rf_PrintValue(sequence);
    NumericVector X = seq2freqProb(sequence);
    // Rf_PrintValue(X);
  //   F<-list()
  //   Q<-list()
  //   QX<-list()
  //   for(o in 1:order) {
  //     F[[o]]<-.seq2matHigh(sequence, o)
  //     Q[[o]]<-sweep(F[[o]], 2, colSums(F[[o]]), '/')
  //     QX[[o]]<-Q[[o]]%*%X
  // # print(QX[[o]])
  //   }
  //   environment(.fn1)=environment()
  //     params<-rep(1/order, order)
  //     model<-Rsolnp::solnp(params, fun = .fn1, eqfun = .eqn1, eqB=1, LB=rep(0, order), control=list(trace=0))
  // # print(model)
  //     lambda=model$pars
  //     return(lambda)
}
