#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(.isProbRcpp)]]
bool isProb(double prob)
{
	if (prob<0 | prob >1) return false;
	return true;
}
