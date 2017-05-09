#include <Rcpp.h>
using namespace Rcpp;
//obtain transition probability matrix from the generator matrix
// [[Rcpp::export]]
NumericMatrix generatorToTransitionMatrix(NumericMatrix gen, bool byrow = true){
  NumericMatrix transMatr(gen.nrow());
  transMatr.attr("dimnames") = gen.attr("dimnames");
  
  if(byrow)
  for(int i = 0; i < gen.nrow(); i++){
    for(int j = 0; j < gen.ncol(); j++){
      if(i != j)
        transMatr(i, j) = -gen(i, j) / gen(i, i);
    }
  }
  else
  for(int j = 0; j < gen.ncol(); j++){
    for(int i = 0; i < gen.nrow(); i++){
      if(i != j)
        transMatr(i, j) = -gen(i, j) / gen(j, j);
    }
  }
  
  return transMatr;
}
