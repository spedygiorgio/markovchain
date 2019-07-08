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
// [[Rcpp::export]]
NumericMatrix generatorToTransitionMatrix(NumericMatrix gen, bool byrow = true){
  NumericMatrix transMatr(gen.nrow());
  transMatr.attr("dimnames") = gen.attr("dimnames");
  
  if (byrow) {
    for (int i = 0; i < gen.nrow(); i++){
      for (int j = 0; j < gen.ncol(); j++){
        if (i != j)
          transMatr(i, j) = -gen(i, j) / gen(i, i);
      }
    }
  } else {
    for (int j = 0; j < gen.ncol(); j++){
      for (int i = 0; i < gen.nrow(); i++){
        if (i != j)
          transMatr(i, j) = -gen(i, j) / gen(j, j);
      }
    }
  }
  
  return transMatr;
}
