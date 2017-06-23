#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <armadillo>
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;
using namespace std;
using namespace RcppArmadillo;
using namespace arma;


// [[Rcpp::export(.ExpectedTimeRCpp)]]
NumericVector ExpectedTimeRcpp(NumericMatrix x,NumericVector y) {
  NumericVector out;
  int size = x.nrow();
  arma::mat T = arma::zeros(size, size);;
  for(int i=0;i<size;i++)
  {
    for(int j=0;j<size;j++)
    {
      T(i,j) = x(i,j);
    }
  }
  arma::vec c = arma::zeros(size);
  for(int i=0;i<size;i++)
    c[i] = y[i];
  out = wrap(solve(T,c));
  return out;
}


