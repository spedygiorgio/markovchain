// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export(.multinomialCIRcpp)]]
List multinomCI(NumericMatrix transMat, NumericMatrix seqMat, double confidencelevel) {
  Function multinomialCI("multinomialCI");
  List res;
  NumericVector v;
  
  // double zscore = stats::qnorm_0(confidencelevel, 1.0, 0.0);
  int nrows = transMat.nrow();
  int ncols = transMat.ncol();
  NumericMatrix lowerEndpointMatr(nrows, ncols);
  NumericMatrix upperEndpointMatr(nrows, ncols);
  double lowerEndpoint, upperEndpoint;
  for(int i = 0; i < nrows; i ++) {
    NumericVector v = seqMat.row(i);
    res = multinomialCI(v, 1-confidencelevel);
    // Rf_PrintValue(res);
    int rsize = res.size();
    for(int j = 0; j < rsize/2; j++) {
      lowerEndpoint = as<double>(res[j]);
      lowerEndpointMatr(i,j%(rsize)) = lowerEndpoint;
      upperEndpoint = as<double>(res[j+rsize/2]);
      upperEndpointMatr(i,j%(rsize)) = upperEndpoint;
    }
  }
  upperEndpointMatr.attr("dimnames") = lowerEndpointMatr.attr("dimnames") = seqMat.attr("dimnames");
  
  List out = List::create(_["confidenceLevel"]=confidencelevel, 
           _["lowerEndpointMatrix"]=lowerEndpointMatr, 
           _["upperEndpointMatrix"]=upperEndpointMatr);
  return out;
}
