// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export(.commclassesKernelRcpp)]]
extern "C" SEXP commclassesKernel(NumericMatrix P){
  unsigned int m = P.ncol(), n;
  CharacterVector stateNames = rownames(P);
  std::vector<int> a;
  arma::vec b, c, d;
  arma::mat T = arma::zeros(m, m);
  unsigned int i = 0;
  int oldSum, newSum, ai;
  while(i < m) {
    a.resize(0);
    a.push_back(i);
    b = arma::zeros<arma::vec>(m);
    b[i] = 1;
    newSum = 0;
    oldSum = 1;
    while(oldSum != newSum) {
      oldSum = 0;
      for(int j = 0; j < b.size(); j ++)
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
      for(int j = 0; j < m; j++) 
        for(int k = 0; k < n; k++)
          c[j] += matr(k, j);
      newSum = 0;
      a.resize(0);
      for(int j = 0; j < b.size(); j++) {
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
  unsigned int sums[tC.n_cols];
  for(unsigned int j = 0; j < T.n_cols; j++) {
    sums[j] = 0;
    for(i = 0; i < T.n_rows; i ++)
      if(tC(i, j) == tT(i, j)) sums[j] ++;
    v[j] = (sums[j] == m);
  }
  C = as<LogicalMatrix>(wrap(Ca));
  C.attr("dimnames") = List::create(stateNames, stateNames);
  v.names() = stateNames;
  return List::create(_["C"] = C, _["v"] = v);
}

//returns the underlying communicating classes
// [[Rcpp::export(.communicatingClassesRcpp)]]
List communicatingClasses(LogicalMatrix adjMatr)
{
  int len = adjMatr.nrow();
  List classesList;
  CharacterVector rnames = rownames(adjMatr);
  for(int i = 0; i < len; i ++) {
    bool isNull = false;
    LogicalVector row2Check = adjMatr(i, _);
    CharacterVector proposedCommClass;
    for(int j = 0; j < row2Check.size(); j++) {
      if(row2Check[j] == true) {
        String rname = rnames[j];
        proposedCommClass.push_back(rname);
      }
    }
    if (i > 0) {
      for(int j = 0; j < classesList.size(); j ++) {
        bool check = false;        
        CharacterVector cv = classesList[j];
        std::set<std::string> s1, s2;
        for(int k = 0; k < cv.size(); k ++) {
          s1.insert(as<std::string>(cv[k]));
          if(proposedCommClass.size() > k) {
            s2.insert(as<std::string>(proposedCommClass[k]));
          }
        }
        check = std::equal(s1.begin(), s1.end(), s2.begin());
        if(check) {
          isNull = true;
          break;
        }
      }
    }
    if(!isNull) 
      classesList.push_back(proposedCommClass);
  }
  return classesList;
}

arma::mat _pow(arma::mat A, int n) {
  arma::mat R = arma::eye(A.n_rows, A.n_rows);
  for(int i = 0; i < n; i ++) 
    R = A*R;
  return R;
}

//communicating states
// [[Rcpp::export(.commStatesFinderRcpp)]]
NumericMatrix commStatesFinder(NumericMatrix matr)
{
  //Reachability matrix
  int dimMatr = matr.nrow();
  arma::mat X(matr.begin(), dimMatr, dimMatr, false);
  arma::mat temp = arma::eye(dimMatr, dimMatr) + arma::sign(X);
  temp = _pow(temp, dimMatr - 1);
  arma::mat m;
  NumericMatrix R = wrap(arma::sign(temp));
  R.attr("dimnames") = List::create(rownames(matr), colnames(matr));
  return R;
}

bool _intersected(CharacterVector v1, CharacterVector v2) {
  CharacterVector::iterator first1 = v1.begin();
  CharacterVector::iterator last1 = v1.end();
  CharacterVector::iterator first2 = v2.begin();
  CharacterVector::iterator last2 = v2.end();
  while(first1!=last1 && first2!=last2) {
    if(*first1 == *first2) return true;
    else if(*first1 < *first2) ++first1;
    else ++first2;    
  }
  return false;
}

// [[Rcpp::export(.summaryKernelRcpp)]]
List summaryKernel(S4 object)
{
  NumericMatrix matr = object.slot("transitionMatrix");
  List temp = commclassesKernel(matr);
  List communicatingClassList = communicatingClasses(temp["C"]);
  List v = temp["v"];
  CharacterVector ns = v.names();
  CharacterVector transientStates; 
  for(int i = 0; i < v.size(); i++) {
    if(v[i] == false)
      transientStates.push_back(ns[i]);
  }
  List closedClasses, transientClasses;

  for(int i = 0; i < communicatingClassList.size(); i ++)
  {
    CharacterVector class2Test = communicatingClassList[i];
    if(_intersected(class2Test,transientStates)) 
        transientClasses.push_back(class2Test); 
      else 
        closedClasses.push_back(class2Test);
  }
  List summaryMc = List::create(_["closedClasses"] = closedClasses,
                                _["transientClasses"] = transientClasses);
  return(summaryMc);
}

//here the kernel function to compute the first passage
// [[Rcpp::export(.firstpassageKernelRcpp)]]
NumericMatrix firstpassageKernel(NumericMatrix P, int i, int n){
  int r = P.nrow(), c = P.ncol();
  arma::mat G(P.begin(), P.nrow(), P.ncol(), false);
  arma::mat Pa = G;
  arma::mat H(n, P.ncol()); //here Thoralf suggestion
  //initializing the first row
  for(int j = 0; j < G.n_cols; j++)
    H(0, j) = G(i-1, j);
  arma::mat E = 1 - arma::eye(P.ncol(), P.ncol());

  for (int m = 1; m < n; m++) {
    G = Pa * (G%E);
    for(int j = 0; j < G.n_cols; j ++) 
      H(m, j) = G(i-1, j);
  }
  NumericMatrix R = wrap(H);
  return R;
}

// greatest common denominator
// [[Rcpp::export(.gcdRcpp)]]
double gcd (int f, int s) {
  int g, n, N, u;
  f = abs(f);
  s = abs(s);
    
  n = std::min(f,s);
  N = std::max(f,s);
  
  if (n==0) {
		g=N;
	}
	else {
		u=1;
		while (u!=0) {
			u=N%n;
			if (u==0) {
				g=n;
			}
			N=n;
			n=u;
		}
	}
	return g;
}
