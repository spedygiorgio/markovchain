#include <Rcpp.h>
using namespace Rcpp;

// From a generator Q build the transition matrix of the EMBEDDED jump chain: from state i the chain
// jumps to j != i with probability -Q[i,j] / Q[i,i]; the diagonal is 0 (the chain always leaves),
// unless the state is absorbing (Q[i,i] = 0), in which case the chain stays with probability 1.
//' @name generatorToTransitionMatrix
//' @title Function to obtain the transition matrix from the generator
//' @description The transition matrix of the embedded DTMC is inferred from the CTMC's generator
//'
//' @usage generatorToTransitionMatrix(gen, byrow = TRUE)
//'
//' @param gen The generator matrix
//' @param byrow Flag to determine if rows (columns) sum to 0
//' @return Returns the transition matrix.
//' 
//' @references
//' Introduction to Stochastic Processes with Applications in the Biosciences (2013), David F.
//' Anderson, University of Wisconsin at Madison
//' 
//' @author Sai Bhargav Yalamanchi
//' @seealso \code{\link{rctmc}},\code{\link{ctmc-class}}
//' @examples
//' energyStates <- c("sigma", "sigma_star")
//' byRow <- TRUE
//' gen <- matrix(data = c(-3, 3, 1, -1), nrow = 2,
//'               byrow = byRow, dimnames = list(energyStates, energyStates))
//' generatorToTransitionMatrix(gen)
//' 
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
