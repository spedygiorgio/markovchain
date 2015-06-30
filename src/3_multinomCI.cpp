// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export(.multinomialCIRcpp)]]
List multinomCI(NumericMatrix transMat, NumericMatrix seqMat, double confidencelevel) {
  Function multinomialCI("multinomialCI");
  List res;
  NumericVector v;
  
  double zscore = stats::qnorm_0(confidencelevel, 1.0, 0.0);
  int nrows = transMat.nrow();
  int ncols = transMat.ncol();
  NumericMatrix lowerEndpointMatr(nrows, ncols);
  NumericMatrix upperEndpointMatr(nrows, ncols);
  double lowerEndpoint, upperEndpoint;
  for(int i = 0; i < nrows; i ++) {
    NumericVector v = seqMat.row(i);
    res = multinomialCI(v, 1-confidencelevel);
//    Rf_PrintValue(res);
    for(int j = 0; j < res.size() - 1; j+=2) {
//      Rf_PrintValue(res[j]);
      lowerEndpoint = as<double>(res[j]);
      lowerEndpointMatr(i,j/2) = lowerEndpoint;
      upperEndpoint = as<double>(res[j+1]);
      upperEndpointMatr(i,j/2) = upperEndpoint;
    }
  }
  upperEndpointMatr.attr("dimnames") = lowerEndpointMatr.attr("dimnames") = seqMat.attr("dimnames");
  
  List out = List::create(_["confidenceLevel"]=confidencelevel, 
           _["lowerEndpointMatrix"]=lowerEndpointMatr, 
           _["upperEndpointMatrix"]=upperEndpointMatr);
  return out;
}

/*** R
library(markovchain)

seq<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
mcfit<-markovchainFit(data=seq,byrow=TRUE)
seqmat<-createSequenceMatrix(seq)
seqmat
.multinomialCIRcpp(mcfit$estimate@transitionMatrix, seqmat, 0.95)
*/
