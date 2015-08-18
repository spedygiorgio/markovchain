#include <Rcpp.h>
using namespace Rcpp;

// double _foo(params) {
//   QX=get("QX")
//   X=get("X")    
//   error=0
//   for (i in 1:length(QX)) {
//     error=error+(params[i]*QX[[i]]-X)
//   }
//   return(sum(error^2))
// }

double _constr(NumericVector params) {
  return(sum(params));
}

//// [[Rcpp::export]]
void fitHigherOrderRcpp(SEXP data, int order = 2) {
  Rcout << "fitHigherOrder " << order << std::endl;
  // Rf_PrintValue(data);
}
// .fitHigherOrder<-function(frequency, sequencelist, order) {
// X=as.numeric(frequency/sum(frequency))
// # ll=vector()
// # Q=alply(.data=seq(1,order,1), .margins=1, .fun=.getQ, sequencelist)
// # transitions=llply(.data=Q, .fun=function(q) q$transition)
// # ll=laply(.data=Q, .fun=function(q) q$ll)
// # QX=llply(.data=transitions, .fun=function(tr) as.matrix(tr)%*%X) 
//   environment(.foo)=environment()
//   params=rep(1/order, order)
//   model=Rsolnp::solnp(pars=params, fun=.foo, eqfun=.constr, eqB=1, LB=rep(0, order), control=list(trace=0))
//   lambda=model$pars
//   return(lambda)
// }

/*** R
*/
