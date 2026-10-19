#include <Rcpp.h>
using namespace Rcpp;

// From a generator Q build the transition matrix of the EMBEDDED jump chain: from state i the chain
// jumps to j != i with probability -Q[i,j] / Q[i,i]; the diagonal is 0 (the chain always leaves),
// unless the state is absorbing (Q[i,i] = 0), in which case the chain stays with probability 1.
//' @name generatorToTransitionMatrix
//' @title Transition matrix of the embedded jump chain
//' @description Computes, from the generator \eqn{Q} of a continuous time
//'   Markov chain, the transition matrix of its embedded (jump) chain: the
//'   probability of jumping from \eqn{i} to \eqn{j \ne i} is
//'   \eqn{-Q_{ij} / Q_{ii}}.
//' @param gen A generator matrix: non-negative off-diagonal rates, with rows
//'   (or columns, if `byrow = FALSE`) summing to zero.
//' @param byrow Logical: are the rates out of each state in its row (default)
//'   or in its column? The result has the same orientation.
//' @return The transition matrix of the embedded chain, with the dimnames of
//'   `gen`. A state with zero exit rate (absorbing) gets probability 1 of
//'   staying where it is.
//' @references Anderson, D. F. (2013). Introduction to Stochastic Processes
//'   with Applications in the Biosciences. University of Wisconsin at
//'   Madison.
//' @author Sai Bhargav Yalamanchi
//' @family continuous time Markov chains
//' @examples
//' energyStates <- c("sigma", "sigma_star")
//' gen <- matrix(c(-3, 3, 1, -1), nrow = 2, byrow = TRUE,
//'               dimnames = list(energyStates, energyStates))
//' generatorToTransitionMatrix(gen)
//' @export
// [[Rcpp::export]]
NumericMatrix generatorToTransitionMatrix(NumericMatrix gen, bool byrow = true){
  // Ho anche corretto la dimensione: usare gen.nrow(), gen.ncol() è più sicuro
  // rispetto a passare un solo argomento al costruttore di Rcpp
  // the loops below read gen(i, i): the generator must be square
  if (gen.nrow() != gen.ncol())
    stop("gen must be a square matrix");
  NumericMatrix transMatr(gen.nrow(), gen.ncol());
  transMatr.attr("dimnames") = gen.attr("dimnames");
  
  // Numbers smaller than tol are treated as zero (absorbing state).
  double tol = 1e-10;
  
  // byrow = TRUE: the rates of leaving i are in row i;  otherwise they are in column i.
  if (byrow) {
    for (int i = 0; i < gen.nrow(); i++){
      if (std::abs(gen(i, i)) < tol) {
        // Issue#211: remaining probability is 1 for absorbing states
        transMatr(i, i) = 1.0; 
      } else {
        for (int j = 0; j < gen.ncol(); j++){
          if (i != j)
            transMatr(i, j) = -gen(i, j) / gen(i, i);
        }
      }
    }
  } else {
    for (int j = 0; j < gen.ncol(); j++){
      if (std::abs(gen(j, j)) < tol) {
        // Issue 211: absorbing state for column-wise logic
        transMatr(j, j) = 1.0; 
      } else {
        for (int i = 0; i < gen.nrow(); i++){
          if (i != j)
            transMatr(i, j) = -gen(i, j) / gen(j, j);
        }
      }
    }
  }
  
  return transMatr;
}
