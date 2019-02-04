// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// returns the column sums of a matrix
NumericVector colSums(NumericMatrix m) {
  NumericVector out;
  for(int i = 0; i < m.cols(); i ++) 
    out.push_back(sum(m.column(i)));
  return out;
}

// poisson distribution
double ppois(double n, double lambda) {
  return R::ppois(n,lambda, true, false);
}

// moments
NumericVector moments(int c, double lambda){
  double a=lambda+c;
  double b=lambda-c;
  double den = 0, poisA = 0, poisB = 0;
  if(b<0) 
    b=0;
  if(b>0)
    den=ppois(a,lambda)-ppois(b-1,lambda);
  if(b==0)
    den=ppois(a,lambda);
  NumericVector mu(4);
  NumericVector mom(5);
  for(int r = 1; r <= 4; r ++){
    poisA=0;
    poisB=0;
    if((a-r) >=0){ poisA=ppois(a,lambda)-ppois(a-r,lambda); }
    if((a-r) < 0){ poisA=ppois(a,lambda); }
    if((b-r-1) >=0){ poisB=ppois(b-1,lambda)-ppois(b-r-1,lambda); }
    if((b-r-1) < 0 && (b-1)>=0){ poisB=ppois(b-1,lambda); }
    if((b-r-1) < 0 && (b-1) < 0){ poisB=0; }
    mu[r - 1]=(pow(lambda,r))*(1-(poisA-poisB)/den);
  }
  mom[0]=mu[0];
  mom[1]=mu[1]+mu[0]-pow(mu[0],2);
  mom[2]=mu[2]+mu[1]*(3-3*mu[0])+(mu[0]-3*pow(mu[0],2)+2*pow(mu[0],3));
  mom[3]=mu[3]+mu[2]*(6-4*mu[0])+mu[1]*(7-12*mu[0]+6*pow(mu[0],2))+mu[0]-4*pow(mu[0],2)+6*pow(mu[0],3)-3*pow(mu[0],4);
  mom[4]=den;
  return mom;
}

// coverage probability for the particular choice of c 
double truncpoi(int c, NumericVector x, double n, int k){
  NumericMatrix m(k,5);
  for(int i = 0; i < k; i ++){
    double lambda=x[i];
    NumericVector mom = moments(c,lambda);
    for(int j = 0; j < 5; j ++)
      m(i,j)=mom[j];
  }
  for(int i = 0; i < k; i ++)
    m(i,3)=m(i,3)-3*m(i,1)*m(i,1);
  
  NumericVector s=colSums(m);
  double s1=s[0];
  double s2=s[1];
  double s3=s[2];
  double s4=s[3];
  
  double probn=1/(ppois(n,n)-ppois(n-1,n));
  double z=(n-s1)/sqrt(s2);
  double g1=s3/(pow(s2,(3.0/2.0)));
  double g2=s4/(pow(s2,2));
  double poly=1.0+g1*(pow(z,3)-3*z)/6.0+g2*(pow(z,4)-6.0*pow(z,2)+3.0)/24.0
              +pow(g1,2)*(pow(z,6)-15.0*pow(z,4)+45.0*pow(z,2)-15.0)/72.0;
  double f=poly*exp(-pow(z,2)/2)/(sqrt(2.0)*R::gammafn(0.5)); 
  double probx=1;
  for(int i = 0; i < k; i++)
    probx=probx*m(i,4);
//   Rcout << poly << " " << g1*(pow(z,3)-3*z)/6+g2*(pow(z,4)-6*pow(z,2)+3)/24 <<
//     " " << 1.0 + g1*(pow(z,3)-3*z)/6+g2*(pow(z,4)-6*pow(z,2)+3)/24 << 
//     " " << pow(g1,2)*(pow(z,6)-15*pow(z,4)+45*pow(z,2)-15)/72 << " " << f << std::endl;
  return(probn*probx*f/sqrt(s2));
}

// multinomial confidence intervals for a row
// [[Rcpp::export(.multinomialCIForRowRcpp)]]
NumericMatrix multinomialCIForRow (NumericVector x, double confidencelevel){
  double n = std::accumulate(x.begin(), x.end(), 0.0);
  int k = x.size();
  double c = 0;
  double p = 0, pold=0;
  for(int cc = 1; cc <= n; cc ++) {
    // Rcout << cc << " " << x << " " << n << " " << k << std::endl;
    p = truncpoi(cc,x,n,k);
    if(p > confidencelevel && pold < confidencelevel) { c = cc; break; };
    pold=p;
  }
  
  NumericMatrix salida(k,2);
  double delta=(confidencelevel-pold)/(p-pold);
  NumericMatrix out(k,5);
  NumericMatrix num(k,1);
  c--;
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
      
    salida(i,0) = out(i,1);
    salida(i,1) = out(i,2);
  }
  return salida;
}

// multinomial confidence intervals
// [[Rcpp::export(.multinomialCIRcpp)]]
List multinomCI(NumericMatrix transMat, NumericMatrix seqMat, double confidencelevel) {
  NumericMatrix res;
  NumericVector v;
  
  int nrows = transMat.nrow();
  int ncols = transMat.ncol();
  NumericMatrix lowerEndpointMatr(nrows, ncols);
  NumericMatrix upperEndpointMatr(nrows, ncols);
  double lowerEndpoint, upperEndpoint;
  for(int i = 0; i < nrows; i ++) {
    v = seqMat.row(i);
    res = multinomialCIForRow(v, confidencelevel);
    for(int j = 0; j < res.rows(); j++) {
      lowerEndpoint = res(j, 0);
      lowerEndpointMatr(i,j) = lowerEndpoint;
      upperEndpoint = res(j, 1);
      upperEndpointMatr(i,j) = upperEndpoint;
    }
  }
  upperEndpointMatr.attr("dimnames") = lowerEndpointMatr.attr("dimnames") = seqMat.attr("dimnames");
  
  List out = List::create(_["confidenceLevel"]=confidencelevel, 
                          _["lowerEndpointMatrix"]=lowerEndpointMatr, 
                          _["upperEndpointMatrix"]=upperEndpointMatr);
  return out;
}
