// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// sequence to frequency probability vector
// [[Rcpp::export]]
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
// [[Rcpp::export]]
NumericMatrix seq2matHigh(CharacterVector sequence, int order) {
  int n = sequence.size();
  CharacterVector states = unique(sequence).sort();
  int nstates = states.length();
  NumericVector colsums(nstates);
  NumericMatrix out(nstates);
  out.attr("dimnames") = List::create(states, states);
  for(int i = 0; i < n - order; i ++) {
    int from = -1, to = -1;
    for (int j = 0; j < nstates; j ++) {
      if(sequence[i] == states[j]) from = j;
      if(sequence[i + order] == states[j]) to = j;
    }
    if(from != -1 && to != -1) {
      out(to, from) ++;
      colsums[from] ++;
    }
  }
  for(int i = 0; i < nstates; i ++) {
    for(int j = 0; j < nstates; j ++)
      out(i, j) /= colsums[j];
  }
  return out;
}

void fn1() {
  // return 1.0;
}

// [[Rcpp::export]]
void fitHigherOrderRcpp(SEXP sequence, int order = 2) {
  NumericVector v = seq2freqProb(sequence);
  arma::vec X(v.begin(), v.size(), false);
  // Rcout << X << std::endl;
  List Q(order), QX(order);
  for(int i = 0; i < order; i ++) {
    NumericMatrix Qi = seq2matHigh(sequence, i + 1);
    arma::mat m(Qi.begin(), Qi.nrow(), Qi.ncol(), false);
    Q[i] = m;
    QX[i] = m * X;
  }
  // Rf_PrintValue(Q);
  // Rf_PrintValue(QX);
  Environment env;
  // env["fn1"] = fn1;
  Rcout << env << std::endl;
  Function solnp("solnp");
  NumericVector params = rep(1.0/order,order);
  // List res = solnp(params, Named("fun", 1));
  //   environment(.fn1)=environment()
  //     params<-rep(1/order, order)
  //     model<-Rsolnp::solnp(params, fun = .fn1, eqfun = .eqn1, eqB=1, LB=rep(0, order), control=list(trace=0))
  // # print(model)
  //     lambda=model$pars
  //     return(lambda)
}
