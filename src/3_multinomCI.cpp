// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export(.multinomialCIRcpp)]]
List multinomCI(NumericMatrix transMat, NumericMatrix seqMat, double confidencelevel) {
  double zscore = stats::qnorm_0(confidencelevel, 1.0, 0.0);
  int nrows = transMat.nrow();
  int ncols = transMat.ncol();
  NumericMatrix lowerEndpointMatr(nrows, ncols);
  NumericMatrix upperEndpointMatr(nrows, ncols);
  // arma::mat tmat(transMat.begin(), transMat.nrow(), transMat.ncol(), false);
  // arma::mat smat(seqMat.begin(), seqMat.nrow(), seqMat.ncol(), false);
  double marginOfError, lowerEndpoint, upperEndpoint;
  for(int i = 0; i < nrows; i ++) {
    for(int j = 0; j < ncols; j ++) {
      marginOfError = zscore * transMat(i, j) / sqrt(seqMat(i, j));
      lowerEndpoint = transMat(i, j) - marginOfError;
      upperEndpoint = transMat(i, j) + marginOfError;
      lowerEndpointMatr(i,j) = (lowerEndpoint > 1.0) ? 1.0 : ((0.0 > lowerEndpoint) ? 0.0 : lowerEndpoint);
      upperEndpointMatr(i,j) = (upperEndpoint > 1.0) ? 1.0 : ((0.0 > upperEndpoint) ? 0.0 : upperEndpoint);
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
