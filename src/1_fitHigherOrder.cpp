#include <Rcpp.h>
using namespace Rcpp;

// sequence to frequency probability vector
NumericVector seq2freqProb (CharacterVector sequence) {
  int n = sequence.size(); 
  CharacterVector states = unique(sequence).sort();
  int nstates = states.length();
  NumericVector v(nstates);
  v.names() = states;
  for(int i = 0; i < n; i ++) {
    v[std::string(sequence[i])] = v[std::string(sequence[i])] + 1.0;
  }
  NumericVector out = v/sum(v);
  out.names() = v.names();
  return out;
}

// sequence to transition matrix for higher order markov chain
NumericMatrix seq2matHigh(CharacterVector sequence, int order) {
  int n = sequence.size();
  CharacterVector states = unique(sequence).sort();
  int nstates = states.length();
  
  NumericMatrix out(nstates);
  out.attr("dimnames") = List::create(states, states);
  for(int i = 0; i < n - order; i ++) {
    int from = -1, to = -1;
    for (int j = 0; j < nstates; j ++) {
      if(sequence[i] == states[j]) from = j;
      if(sequence[i + order] == states[j]) to = j;
    }
    if(from != -1 && to != -1)
      out(to, from) ++;
    // out(to, from) = out(to, from) + 1.0;
  }
  // Rf_PrintValue(out);
  return out;
}
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
  List F(order), Q(order), QX(order);
  for(int i = 0; i < order; i ++) {
    F[i] = seq2matHigh(sequence, i + 1);
  }
  // Rf_PrintValue(F);
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
