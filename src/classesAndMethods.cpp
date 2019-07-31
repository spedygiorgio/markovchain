// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// check if prob is probability or not
// [[Rcpp::export(.isProbability)]]
bool isProb(double prob) {
  return (prob >= 0 && prob <= 1);
}

// doubt
// [[Rcpp::export(.isGenRcpp)]]
bool isGen(NumericMatrix gen) {
  for (int i = 0; i < gen.nrow(); i++)
    for (int j = 0; j < gen.ncol(); j++)
      if ((i == j && gen(i, j) > 0) || (i != j && gen(i, j) < 0))  
        return false;

  return true;
}

SEXP commClassesKernel(NumericMatrix P);

// method to convert into canonic form a markovchain object
// [[Rcpp::export(.canonicFormRcpp)]]
SEXP canonicForm(S4 object) {
  NumericMatrix P = object.slot("transitionMatrix");
  List comclasList = commClassesKernel(P);
  LogicalVector vu = comclasList["closed"];
  NumericVector u, w; 

  for (int i = 0; i < vu.size(); i ++) {
    if(vu[i]) u.push_back(i);
    else w.push_back(i);
  }
  
  LogicalMatrix Cmatr = comclasList["classes"];
  NumericVector R, p;
  LogicalVector crow;
  
  while (u.size() > 0) {
    R.push_back(u[0]);
    crow = Cmatr(u[0], _);

    for (int i = 0; i < crow.size(); i++) 
      vu[i] = vu[i] * !crow[i];

    u = NumericVector::create();

    for (int i = 0; i < vu.size(); i ++) 
      if(vu[i]) u.push_back(i);
  }
  
  for (int i = 0; i < R.size(); i ++) {
    crow = Cmatr(R[i], _);

    for (int j = 0; j < crow.size(); j++) 
      if(crow[j]) p.push_back(j);
  }
  
  for (NumericVector::iterator it = w.begin(); it != w.end(); it++)
    p.push_back(*it);
  
  NumericMatrix Q(p.size());
  CharacterVector rnames(P.nrow());
  CharacterVector cnames(P.ncol());
  CharacterVector r = rownames(P);
  CharacterVector c = colnames(P);

  for (int i = 0; i < p.size(); i ++) {
    rnames[i] = r[p[i]];
    for (int j = 0; j < p.size(); j ++) {
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
  
  for (int i=0; i<m.nrow(); i++)
    for (int j=0; j<m.ncol(); j++)
      x[i][j] = m(i,j);
  
  sort(x.begin(), x.end());
  
  return(wrap(x));
}

inline bool approxEqual(const double& a, const double& b) {
  if (a >= b)
    return (a - b) <= 1E-7;
  else
    return approxEqual(b, a);
}

inline bool approxEqual(const cx_double& a, const cx_double& b){
  double x = a.real() - b.real();
  double y = a.imag() - b.imag();
  
  return (x*x - y*y) <= 1E-14;
}

mat mcEigen(NumericMatrix t) {
  cx_mat transitionMatrix = as<cx_mat>(t);
  cx_vec eigvals;
  cx_mat eigvecs;
  // 1 + 0i
  cx_double cx_one(1.0, 0);

  // If transition matrix is hermitian (symmetric in real case), use
  // more efficient implementation to get the eigenvalues and vectors
  if (transitionMatrix.is_hermitian()) {
    vec real_eigvals;
    eig_sym(real_eigvals, eigvecs, transitionMatrix);
    eigvals.resize(real_eigvals.size());
    
    // eigen values are real, but we need to cast them to complex values
    // to perform the rest of the algorithm
    for (int i = 0; i < real_eigvals.size(); ++i)
      eigvals[i] = cx_double(real_eigvals[i], 0);
  } else
    eig_gen(eigvals, eigvecs, transitionMatrix, "balance");

  std::vector<int> whichOnes;
  std::vector<double> colSums;
  double colSum;
  mat real_eigvecs = real(eigvecs);
  int numRows = real_eigvecs.n_rows;
    
  // Search for the eigenvalues which are 1 and store 
  // the sum of the corresponding eigenvector
  for (int j = 0; j < eigvals.size(); ++j) {
    if (approxEqual(eigvals[j], cx_one)) {
      whichOnes.push_back(j);
      colSum = 0;
      
      for (int i = 0; i < numRows; ++i)
        colSum += real_eigvecs(i, j);
      
      colSums.push_back(colSum);
    }
  }

  // Normalize eigen vectors
  int numCols = whichOnes.size();
  mat result(numRows, numCols);
  
  for (int i = 0; i < numRows; ++i)
    for (int j = 0; j < numCols; ++j)
      result(i, j) = real_eigvecs(i, whichOnes[j]) / colSums[j];


  return result;
}

