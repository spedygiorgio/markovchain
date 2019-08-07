// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// sequence to frequency probability vector
//' @export
// [[Rcpp::export]]
NumericVector seq2freqProb(CharacterVector sequence) {
  int n = sequence.size(); 
  CharacterVector states = unique(sequence).sort();
  int nstates = states.length();
  NumericVector v(nstates);
  v.names() = states;
  
  for (int i = 0; i < n; i ++)
    v[std::string(sequence[i])] = v[std::string(sequence[i])] + 1.0;
  
  NumericVector out = v/sum(v);
  out.names() = v.names();
  
  return out;
}

// sequence to transition matrix for higher order markov chai
//' @export
// [[Rcpp::export]]
NumericMatrix seq2matHigh(CharacterVector sequence, int order) {
  int n = sequence.size();
  CharacterVector states = unique(sequence).sort();
  int nstates = states.length();
  NumericVector colsums(nstates);
  NumericMatrix out(nstates);
  out.attr("dimnames") = List::create(states, states);
  
  for (int i = 0; i < n - order; i++) {
    int from = -1, to = -1;

    for (int j = 0; j < nstates; j++) {
      if (sequence[i] == states[j])
        from = j;
      if (sequence[i + order] == states[j])
        to = j;
    }
    
    if (from != -1 && to != -1) {
      out(to, from) ++;
      colsums[from] ++;
    }
  }
  
  for (int i = 0; i < nstates; i ++) {
    for (int j = 0; j < nstates; j ++)
      out(i, j) /= colsums[j];
  }
  
  return out;
}
