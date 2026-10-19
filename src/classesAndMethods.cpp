// This file is compiled C++ code called from R through Rcpp. Lines such as
// // [[Rcpp::export(.name)]] make Rcpp generate an R function .name that calls the C++ function
// (the glue code is in RcppExports.cpp / R/RcppExports.R, generated automatically: do not edit it).
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <functional>
#include <unordered_map>
#include <string>
using namespace Rcpp;
using namespace arma;
using namespace std;


// Checks that a matrix is a valid GENERATOR of a continuous time chain: diagonal entries <= 0 and
// off-diagonal entries >= 0 (rows sum to zero is checked in R).
// TODO meaning of this method
// [[Rcpp::export(.isGenRcpp)]]
bool isGen(NumericMatrix gen) {
  for (int i = 0; i < gen.nrow(); i++)
    for (int j = 0; j < gen.ncol(); j++)
      if ((i == j && gen(i, j) > 0) || (i != j && gen(i, j) < 0))  
        return false;

  return true;
}
