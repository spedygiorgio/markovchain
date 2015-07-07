// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// double moments(c,lambda){
// a=lambda+c;
// b=lambda-c;
// if(b<0){ 
// b=0;
// }
// if(b>0){ 
// # den=poisson(lambda,a)-poisson(lambda,b-1);
// den=ppois(a,lambda)-ppois(b-1,lambda);
// }
// if(b==0){
// # den=poisson(lambda,a);
// den=ppois(a,lambda);
// }
// mu=mat.or.vec(4,1); 
// # mom es global y se usa fuera de esta funci?n
// mom=mat.or.vec(5,1); 
// for(r in 1:4){
// poisA=0;
// poisB=0;
// if((a-r) >=0){ poisA=ppois(a,lambda)-ppois(a-r,lambda); }
// if((a-r) < 0){ poisA=ppois(a,lambda); }
// if((b-r-1) >=0){ poisB=ppois(b-1,lambda)-ppois(b-r-1,lambda); }
// if((b-r-1) < 0 && (b-1)>=0){ poisB=ppois(b-1,lambda); }
// if((b-r-1) < 0 && (b-1) < 0){ poisB=0; }
// mu[r]=(lambda^r)*(1-(poisA-poisB)/den);
// }
// mom[1]=mu[1];
// mom[2]=mu[2]+mu[1]-mu[1]^2;
// mom[3]=mu[3]+mu[2]*(3-3*mu[1])+(mu[1]-3*mu[1]^2+2*mu[1]^3);
// mom[4]=mu[4]+mu[3]*(6-4*mu[1])+mu[2]*(7-12*mu[1]+6*mu[1]^2)+mu[1]-4*mu[1]^2+6*mu[1]^3-3*mu[1]^4;
// mom[5]=den;
// mom
// }

double truncpoi(int c, NumericVector x, double n, int k){
  NumericMatrix m(k,5);
  for(int i = 0; i < k; i ++){
    double lambda=x[i];
//     NumericVector mom = moments(c,lambda);
//     for(int j = 0; j < 5; j ++)
//       m[i,j]=mom[j];
  }
  for(int i = 0; i < k; i ++)
    m(i,3)=m(i,3)-3*m(i,1)*m(i,1);
    
//   NumericVector s=colSums(m);
//   double s1=s[0];
//   double s2=s[1];
//   double s3=s[2];
//   double s4=s[3];
  //   
  //   probn=1/(ppois(n,n)-ppois(n-1,n));
  //   z=(n-s1)/sqrt(s2);
  //   g1=s3/(s2^(3/2));
  //   g2=s4/(s2^2);
  //   poly=1+g1*(z^3-3*z)/6+g2*(z^4-6*z^2+3)/24
  //   +g1^2*(z^6-15*z^4+45*z^2-15)/72;
  //   f=poly*exp(-z^2/2)/(sqrt(2)*gamma(0.5));
  //   probx=1;
  //   for(i in 1:k){
  //     probx=probx*m[i,5];
  //   }
  //   return(probn*probx*f/sqrt(s2));
  return 0;
}

NumericMatrix multinomialCIForRow (NumericVector x, double confidencelevel){
  double n = std::accumulate(x.begin(), x.end(), 0.0);
  int k = x.size();
  // Rf_PrintValue(x);
  // NumericVector pp = x/n;
  double c = 0;
  double p = 0, pold=0;
  for(int cc = 1; cc <= n; cc ++) {
    p = truncpoi(cc,x,n,k);
    if(p > confidencelevel && pold < confidencelevel) { c = cc; break; };
    pold=p;
  }
  
  NumericMatrix salida(k,2);
  double delta=(confidencelevel-pold)/(p-pold);
  NumericMatrix out(k,5);
  NumericMatrix num(k,1);
  c--;
  double vol1=1;
  double vol2=1;
  double obsp = 0;
  for(int i = 0; i < k; i++) {
    num(i, 0) = i;
    obsp=x[i]/n;
    out(i,0)=obsp;
    out(i,1)=obsp-c/n;
    out(i,2)=obsp+c/n+2*delta/n;
    if(out(i,1)<0){ out(i,1)=0; }
    if(out(i,2)>1){ out(i,2)=1; }
    out(i,3)=obsp-c/n-1/n;
    out(i,4)=obsp+c/n+1/n;
//     if(out(i,1)<0){ out(i,1)=0; }
//     if(out(i,2)>1){ out(i,2)=1; }
//     vol1=vol1*(out(i,2)-out(i,1));
//     vol2=vol2*(out(i,4)-out(i,3));
      
    salida(i,0) = out(i,1);
    salida(i,1) = out(i,2);
  }
  //   c1=c('PROPORTION', 'LOWER(SG)','UPPER(SG)','LOWER(C+1)','UPPER(C+1)');
  //   cov=100*(1-alpha);
  //   sg=(x+delta)/n;
  //   c2=c('SG-midpoint');
  return salida;
}

// [[Rcpp::export(.multinomialCIRcpp2)]]
List multinomCI2(NumericMatrix transMat, NumericMatrix seqMat, double confidencelevel) {
  NumericMatrix res;
  NumericVector v;
  
  int nrows = transMat.nrow();
  int ncols = transMat.ncol();
  NumericMatrix lowerEndpointMatr(nrows, ncols);
  NumericMatrix upperEndpointMatr(nrows, ncols);
  double lowerEndpoint, upperEndpoint;
  for(int i = 0; i < nrows; i ++) {
    v = seqMat.row(i);
    res = multinomialCIForRow(v, 1-confidencelevel);
    // Rf_PrintValue(res);
    // int rsize = res.size();
    //     for(int j = 0; j < rsize/2; j++) {
    //       lowerEndpoint = as<double>(res[j]);
    //       lowerEndpointMatr(i,j%(rsize)) = lowerEndpoint;
    //       upperEndpoint = as<double>(res[j+rsize/2]);
    //       upperEndpointMatr(i,j%(rsize)) = upperEndpoint;
    //     }
  }
  upperEndpointMatr.attr("dimnames") = lowerEndpointMatr.attr("dimnames") = seqMat.attr("dimnames");
  
  List out = List::create(_["confidenceLevel"]=confidencelevel, 
                          _["lowerEndpointMatrix"]=lowerEndpointMatr, 
                          _["upperEndpointMatrix"]=upperEndpointMatr);
  return out;
}

// [[Rcpp::export(.multinomialCIRcpp)]]
List multinomCI(NumericMatrix transMat, NumericMatrix seqMat, double confidencelevel) {
  Function multinomialCI("multinomialCI");
  List res;
  NumericVector v;
  
  int nrows = transMat.nrow();
  int ncols = transMat.ncol();
  NumericMatrix lowerEndpointMatr(nrows, ncols);
  NumericMatrix upperEndpointMatr(nrows, ncols);
  double lowerEndpoint, upperEndpoint;
  for(int i = 0; i < nrows; i ++) {
    v = seqMat.row(i);
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
