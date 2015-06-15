// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export(.isProbRcpp)]]
bool isProb(double prob)
{
	if (prob<0 || prob >1) return false;
	return true;
}

extern "C" SEXP commclassesKernel2(NumericMatrix P){
  unsigned int m = P.ncol(), n;
  CharacterVector stateNames = rownames(P);
  std::vector<int> a;
  arma::vec b, c, d;
  arma::mat T = arma::zeros(m, m);
  unsigned int i = 0;
  int oldSum, newSum; // SAI fix
  while(i < m) {
    a.resize(0);
    a.push_back(i);
    b = arma::zeros<arma::vec>(m);
    b[i] = 1;
    newSum = 0;
    oldSum = 1;
    while(oldSum != newSum) {
      oldSum = 0;
      for(unsigned int j = 0; j < b.size(); j ++)
        if(b[j] > 0) oldSum += (j + 1);
      n = a.size();
      NumericVector temp; 
      NumericMatrix matr(n, m);
      for(unsigned int j = 0; j < n; j ++) {
        temp = P.row(a[j]);
        for(int k = 0; k < temp.size(); k++) 
          matr(j, k) = temp[k];
      }
      c = arma::zeros<arma::vec>(m);
      for(unsigned int j = 0; j < m; j++) 
        for(unsigned int k = 0; k < n; k++)
          c[j] += matr(k, j);
      newSum = 0;
      a.resize(0);
      for(unsigned int j = 0; j < b.size(); j++) {
        if(c[j] > 0) {
          b[j] = 1; a.push_back(j);
        }
        if(b[j] > 0) newSum += (j + 1);
      }
    }
    for(unsigned int j = 0; j < b.size(); j ++)
      T(i, j) = b[j];
    i++;
  }
  arma::mat F = arma::trans(T);
  LogicalMatrix C;
  arma::mat Ca(T.n_rows, T.n_cols);
  for(i = 0; i < T.n_rows; i ++) {
   for(unsigned int j = 0; j < T.n_cols; j++) {
      Ca(i, j) = (T(i, j) > 0 && F(i, j) > 0);
   }
  }
  LogicalVector v(T.n_cols);
  arma::mat tC = Ca.t();
  arma::mat tT = T.t();
  IntegerVector sums(tC.n_cols);
  for(unsigned int j = 0; j < T.n_cols; j++) {
    sums[j] = 0;
    for(i = 0; i < T.n_rows; i ++)
      if(tC(i, j) == tT(i, j)) sums[j] ++;
    v[j] = (sums[j] == (int)m);
  }
  C = as<LogicalMatrix>(wrap(Ca));
  C.attr("dimnames") = List::create(stateNames, stateNames);
  v.names() = stateNames;
  return List::create(_["C"] = C, _["v"] = v);
}

// method to convert into canonic form a markovchain object
// [[Rcpp::export(.canonicFormRcpp)]]
extern "C" SEXP canonicForm (S4 object)
{
  NumericMatrix P = object.slot("transitionMatrix");
  List comclasList = commclassesKernel2(P);
  LogicalVector vu = comclasList["v"];
  NumericVector u, w; 
  for(int i = 0; i < vu.size(); i ++) {
    if(vu[i]) u.push_back(i);
    else w.push_back(i);
  }
  
  LogicalMatrix Cmatr = comclasList["C"];
  NumericVector R, p;
  LogicalVector crow;
  while(u.size()>0)
  {
    R.push_back(u[0]);
    crow = Cmatr(u[0], _);
    for(int i = 0; i < crow.size(); i++) 
      vu[i] = vu[i] * !crow[i];
    u = NumericVector::create();
    for(int i = 0; i < vu.size(); i ++) 
      if(vu[i]) u.push_back(i);
  }
  for (int i = 0; i < R.size(); i ++)
  {
    crow = Cmatr(R[i], _);
    for(int j = 0; j < crow.size(); j++) 
      if(crow[j]) p.push_back(j);
  }
  for(NumericVector::iterator it = w.begin(); it != w.end(); it++)
    p.push_back(*it);
  NumericMatrix Q(p.size());
  CharacterVector rnames(P.nrow());
  CharacterVector cnames(P.ncol());
  CharacterVector r = rownames(P);
  CharacterVector c = colnames(P);
  for(int i = 0; i < p.size(); i ++) {
    rnames[i] = r[p[i]];
    for(int j = 0; j < p.size(); j ++) {
      Q(i, j) = P(p[i], p[j]);
      cnames[j] = c[p[j]];
    }
  }
  Q.attr("dimnames") = List::create(rnames, cnames);
  S4 out("markovchain"); 
  out.slot("transitionMatrix") = Q;
  out.slot("name") = object.slot("name");
  return out;
}

