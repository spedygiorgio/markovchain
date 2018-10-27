// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// check if prob is probability or not
// [[Rcpp::export(.isProbability)]]
bool isProb(double prob) {
  prob >= 0 && prob <= 1;
}
// doubt
// [[Rcpp::export(.isGenRcpp)]]
bool isGen(NumericMatrix gen) {
  for(int i = 0; i < gen.nrow(); i++)
    for(int j = 0; j < gen.ncol(); j++)
      if((i == j && gen(i, j) > 0) || (i != j && gen(i, j) < 0))  
        return false;
  return true;
}

SEXP commclassesKernel(NumericMatrix P);

// method to convert into canonic form a markovchain object
// [[Rcpp::export(.canonicFormRcpp)]]
SEXP canonicForm (S4 object)
{
  NumericMatrix P = object.slot("transitionMatrix");
  List comclasList = commclassesKernel(P);
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




// Function to sort a list of vectors lexicographically
// Input should be given as a matrix where each row represents a vector
// [[Rcpp::export(.lexicographical_sort)]]
SEXP lexicographicalSort(SEXP y) {
  NumericMatrix m(y);
  std::vector< std::vector<double> > x(m.nrow(), std::vector<double>(m.ncol()));
  
  for(int i=0; i<m.nrow(); i++)
    for(int j=0; j<m.ncol(); j++)
      x[i][j] = m(i,j);
  
  sort(x.begin(), x.end());
  
  return(wrap(x));
}
