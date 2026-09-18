// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// sequence to frequency probability vector
//' @export
// [[Rcpp::export]]
NumericVector seq2freqProb(CharacterVector sequence) {
  if (sequence.size() < 1) stop("sequence must not be empty");
  for (R_xlen_t i = 0; i < sequence.size(); ++i)
    if (CharacterVector::is_na(sequence[i]))
      stop("sequence must not contain missing values");
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

// sequence to transition matrix for higher order markov chain.
//
// SUMMARY FOR REVIEWERS: this builds the lag-`order` marginal transition
// matrix out(to, from) = P(state at t+order = to | state at t = from),
// used by fitHigherOrder() (see fitHigherOrder.R) to implement the Ching,
// Ng & Fung (2008) higher-order Markov chain model: a weighted combination
// of these lag-k matrices for k = 1..order, with weights fit by quadratic
// programming (Rsolnp::solnp) to best reproduce the overall stationary
// distribution. Each lag matrix is column-stochastic (columns are the
// "from" state, matching how fitHigherOrder.R multiplies Q[[o]] %*% X).
//
// FIX (see git history / PR notes for details): previously, a state that
// never occurred as a "from" state at a given lag (colsums[j] == 0, e.g.
// a rare state, or a lag close to the sequence length) produced an entire
// column of NaN (0/0) with no warning. Because Q %*% X mixes all rows, a
// single NaN column silently turned the ENTIRE product -- and hence the
// entire fitHigherOrder() objective function -- into NaN, making the fit
// fail with no indication why. The fallback below (uniform distribution
// for that column) fixes this and matches the "no data available" 
// convention used elsewhere in the package (e.g. generateCI's handling of
// an empty row in fittingFunctions.cpp).
//' @export
// [[Rcpp::export]]
NumericMatrix seq2matHigh(CharacterVector sequence, int order) {
  const int n = sequence.size();
  if (n < 2) stop("sequence must contain at least two observations");
  for (R_xlen_t i = 0; i < sequence.size(); ++i)
    if (CharacterVector::is_na(sequence[i]))
      stop("sequence must not contain missing values");
  if (order < 1 || order >= n)
    stop("order must be positive and smaller than the sequence length");
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
  
  // A column with colsums[j] == 0 means state j never occurred as the
  // "from" state at this lag (e.g. a rare state, or a lag close to the
  // sequence length) -- dividing by zero would silently fill that whole
  // column with NaN, which then contaminates every entry of Q %*% X in
  // fitHigherOrder() (a single NaN column turns the entire product into
  // NaN), causing the quadratic-programming fit to fail with no clear
  // indication why. Falling back to a uniform distribution for that
  // column instead matches the convention used elsewhere in the package
  // for states with no data (see e.g. the `sanitize` handling in
  // createSequenceMatrix / generateCI).
  for (int i = 0; i < nstates; i ++) {
    for (int j = 0; j < nstates; j ++) {
      if (colsums[j] == 0)
        out(i, j) = 1.0 / nstates;
      else
        out(i, j) /= colsums[j];
    }
  }
  
  return out;
}
