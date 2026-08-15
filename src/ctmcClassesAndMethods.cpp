#include <Rcpp.h>
using namespace Rcpp;

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
// [[Rcpp::export(.generatorToTransitionMatrix)]]
NumericMatrix generatorToTransitionMatrix(NumericMatrix gen, bool byrow = true){
  // Ho anche corretto la dimensione: usare gen.nrow(), gen.ncol() è più sicuro
  // rispetto a passare un solo argomento al costruttore di Rcpp
  NumericMatrix transMatr(gen.nrow(), gen.ncol());
  transMatr.attr("dimnames") = gen.attr("dimnames");
  
  double tol = 1e-10;
  
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
