// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <functional>
#include <unordered_map>
#include <string>
using namespace Rcpp;
using namespace arma;
using namespace std;


// TODO meaning of this method
// [[Rcpp::export(.isGenRcpp)]]
bool isGen(NumericMatrix gen) {
  for (int i = 0; i < gen.nrow(); i++)
    for (int j = 0; j < gen.ncol(); j++)
      if ((i == j && gen(i, j) > 0) || (i != j && gen(i, j) < 0))  
        return false;

  return true;
}