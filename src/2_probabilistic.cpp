// [[Rcpp::depends(RcppArmadillo)]]

//#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export(.commclassesKernelRcpp)]]
extern "C" SEXP commclassesKernel(NumericMatrix P){
  unsigned int m = P.ncol(), n;
  CharacterVector stateNames = rownames(P);
  arma::vec a, b, c, d;
  arma::mat T = arma::zeros(m);
  unsigned int oldSum, newSum, i = 0;
  while(i <= m) {
    a = i;
    b = arma::vec(m);
    b[i] = 1;
    newSum = 0;
    oldSum = 1;
    while(oldSum != newSum) {;
      oldSum = sum(find(b > 0));
      n = a.size();
      NumericVector r, temp; 
      for(unsigned int j = 0; j < n; j ++) {
        temp = P.row(a[j]);
        for(NumericVector::iterator it = temp.begin(); it != temp.end(); it ++)
          r.push_back(*it);
      }
      NumericMatrix matr(n, m, r.begin());
      c = sum(matr);
      arma::vec d = c;
  		n = d.size();
      for(arma::vec::iterator it = d.begin(); it != d.end(); it++)
        b[*it] = 1;
      newSum = 0;
      for(arma::vec::iterator it = b.begin(); it != b.end(); it ++)
        if(*it > 0) newSum += *it;
      newSum = sum(find(b > 0));
	    a = d;
    }
    T.insert_rows(i, b);
    i++;
  }

  arma::mat F = arma::trans(T);
  NumericMatrix C;
  arma::mat Ca(T.n_rows, T.n_cols);
  for(i = 0; i < T.n_rows; i ++)
    for(unsigned int j = 0; j < T.n_cols; j++)
      Ca(i, j) = (T(i, j) > 0 && F(i, j) > 0);
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
  for(int i = 0; i < len; i ++) {
    LogicalVector row2Check = adjMatr(i, _);
    CharacterVector rnames = row2Check.names();
    CharacterVector proposedCommClass;// = names(which(row2Check == true));
    for(int j = 0; j < row2Check.size(); j++) 
      if(row2Check[j] == true) 
          proposedCommClass.push_back(rnames(j));
    if (i > 0) {
      for(int j = 0; j < classesList.size(); j ++) {
        bool check = false;        
        CharacterVector cv = classesList[j];
        std::set<std::string> s1, s2;
        for(int k = 0; k < cv.size(); k ++) {
          s1.insert(as<std::string>(cv[k]));
          s2.insert(as<std::string>(proposedCommClass[k]));
        }
        check = std::equal(s1.begin(), s1.end(), s2.begin());
        if(check) {
          proposedCommClass = R_NilValue; break;
        }
      }
    }
    if(!Rf_isNull(proposedCommClass) && proposedCommClass.size() > 0) 
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
  arma::mat R = arma::sign(temp);
  return wrap(R);
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
  CharacterVector transientStates; //<-names(which(temp$v==FALSE))
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
  arma::mat G(P.begin(), P.nrow(), P.ncol());
  arma::mat Pa = G;
  arma::mat H(n, P.ncol()); //here Thoralf suggestion
  H.insert_rows(0, G.row(0)); //initializing the first row  
  arma::mat E = 1 - arma::eye(P.ncol(), P.ncol());
  
  for (int m = 1; m < n; m++) {
    G = Pa * (G*E);
    //H<-rbind(H,G[i,]) //removed thanks to Thoralf 
    H.insert_rows(m, G.row(i)); //here Thoralf suggestion
  }
  return wrap(H);
}

// greatest common denominator: to be moved in Rcpp
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

