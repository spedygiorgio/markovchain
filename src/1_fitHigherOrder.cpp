#include <Rcpp.h>
using namespace Rcpp;

NumericVector _rotate(NumericVector x, int n) {
  int l=x.size();
  n=n % l;
  // n=n %% l;
  if (n == 0) {
    return(x);
  }
//   tmp=x[(l-n+1):l]
//   x[(n+1):l]=x[1:(l-n)]
//   x[1:n]=tmp
    return(x);
}

// List _getQ=function(i, sequenceList) { 
// clicks=sequenceList
//   for (j in 1:i) {
//     clicks=rbind(clicks, "[[-]]")
//   }
//   clicks=unlist(clicks, use.names=F)
//     clicks2=.rotate(clicks, -i)
// # dat=data.table(clicks, clicks2) 
// #   transition=as.data.frame(dcast.data.table(dat, clicks~clicks2, fun.aggregate=length, value.var="clicks2"))
// #   transition=transition[,-1]  
// #   pos=which(names(transition)=="[[-]]")
// #   rnames=names(transition)[-pos]
// #   transition=transition[,-pos]
// #   transition=transition[-pos,]
// #   sums=colSums(t(transition))
// #   sums[sums==0]=1
// #   ll=sum(transition*log(transition/sums), na.rm=T)
// #   transition=as.data.frame(t(transition/sums))
// #   names(transition)=rnames
// #   rownames(transition)=rnames
// #   return(list(ll=ll, transition=transition))
// }

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

// [[Rcpp::export]]
void fitHigherOrder(SEXP data, int order = 2) {
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
