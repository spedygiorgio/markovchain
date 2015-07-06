// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;

template <typename T>
T sortByDimNames(const T m);

// [[Rcpp::export(.commclassesKernelRcpp)]]
SEXP commclassesKernel(NumericMatrix P){
  unsigned int m = P.ncol(), n;
  CharacterVector stateNames = rownames(P);
  std::vector<int> a;
  arma::vec b, c, d;
  arma::mat T = arma::zeros(m, m);
  unsigned int i = 0;
  int oldSum, newSum;
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
      c = arma::zeros<arma::vec>(m);
      for(unsigned int j = 0; j < n; j ++) {
        temp = P.row(a[j]);
        for(int k = 0; k < temp.size(); k++)  {
          matr(j, k) = temp[k];
          c[k] += matr(j, k);
        }
      }
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

//returns the underlying communicating classes
// [[Rcpp::export(.communicatingClassesRcpp)]]
List communicatingClasses(S4 object)
{
  NumericMatrix matr = object.slot("transitionMatrix");
  List temp = commclassesKernel(matr);
  LogicalMatrix adjMatr = temp["C"];
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

// [[Rcpp::export(.recurrentClassesRcpp)]]
List recurrentClasses(S4 object)
{
  NumericMatrix matr = object.slot("transitionMatrix");
  List temp = commclassesKernel(matr);
  LogicalMatrix adjMatr = temp["C"];
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
      } else if(matr(i, j) > 0) {
        isNull = true;
        break;
      }
    }
    if (i > 0) {
      for(int j = 0; j < classesList.size(); j ++) {
        bool check = false;        
        CharacterVector cv = classesList[j];
        std::set<std::string> s1, s2;
        for(int k = 0; k < cv.size(); k ++) {
          s1.insert(as<std::string>(cv[k]));
          if(proposedCommClass.size() > k) 
            s2.insert(as<std::string>(proposedCommClass[k]));
        }
        if(!s1.empty() && !s2.empty())
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
  NumericMatrix R = wrap(arma::sign(temp));
  R.attr("dimnames") = matr.attr("dimnames");
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
  List communicatingClassList = communicatingClasses(object);
  // List communicatingClassList = communicatingClasses(temp["C"]);
  List v = temp["v"];
  CharacterVector ns = v.names();
  CharacterVector transientStates; 
  for(int i = 0; i < v.size(); i++) {
    if(bool(v[i]) == false)
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
  arma::mat G = as<arma::mat>(P);
  arma::mat Pa = G;
  arma::mat H(n, P.ncol()); //here Thoralf suggestion
  //initializing the first row
  for(unsigned int j = 0; j < G.n_cols; j++)
    H(0, j) = G(i-1, j);
  arma::mat E = 1 - arma::eye(P.ncol(), P.ncol());

  for (int m = 1; m < n; m++) {
    G = Pa * (G%E);
    for(unsigned int j = 0; j < G.n_cols; j ++) 
      H(m, j) = G(i-1, j);
  }
  NumericMatrix R = wrap(H);
  return R;
}

// greatest common denominator
// [[Rcpp::export(.gcdRcpp)]]
int gcd (int a, int b) {
  int c;
  a = abs(a);
  b = abs(b);
  while ( a != 0 ) {
    c = a; a = b%a;  b = c;
  }
  return b;
}

//function to  get the period of a DTMC
// [[Rcpp::export(period)]]
int period(S4 object) {
  Function isIrreducible("is.irreducible");
  List res = isIrreducible(object);
  if(!res[0]) {
    warning("The matrix is not irreducible");
    return 0;
  } else {
    NumericMatrix P = object.slot("transitionMatrix");
    int n = P.ncol();
    arma::vec v(n);
    std::vector<double> r, T(1), w;
    v[0] = 1;
    int d = 0, m = T.size(), i = 0, j = 0;
    while(m>0 && d!=1) {
      i = T[0];
      T.erase(T.begin());
      w.push_back(i);
      j = 0;
      while(j < n) {
        if(P(i,j) > 0) {
          r.insert(r.end(), w.begin(), w.end());
          r.insert(r.end(), T.begin(), T.end());
          double k = 0;
          for(std::vector<double>::iterator it = r.begin(); it != r.end(); it ++) 
            if(*it == j) k ++;
          if(k > 0) {
             int b = v[i] + 1 - v[j];
             d = gcd(d, b);
          } else {
            T.push_back(j);
            v[j] = v[i] + 1;
          }
        }
        j ++;
      }
      m = T.size();
    }
    v = v - floor(v/d)*d;
    return d;
  }
}

// [[Rcpp::export]]
double predictiveDistribution(CharacterVector stringchar, CharacterVector newData, NumericMatrix hyperparam = NumericMatrix()) {
  // construct list of states
  CharacterVector elements = stringchar;
  for(int i = 0; i < newData.size(); i++)
    elements.push_back(newData[i]);
  
  elements = unique(elements).sort();
  int sizeMatr = elements.size();
  
  // if no hyperparam argument provided, use default value of 1 for all 
  if(hyperparam.nrow() == 1 && hyperparam.ncol() == 1){
    NumericMatrix temp(sizeMatr, sizeMatr);
    temp.attr("dimnames") = List::create(elements, elements);
    for(int i = 0; i < sizeMatr; i++)
      for(int j = 0; j < sizeMatr; j++)
        temp(i, j) = 1;
    hyperparam = temp;
  }
  
  // validity check
  if(hyperparam.nrow() != hyperparam.ncol())
    stop("Dimensions of the hyperparameter matrix are inconsistent");
    
  if(hyperparam.nrow() < sizeMatr)
    stop("Hyperparameters for all state transitions must be provided");
    
  List dimNames = hyperparam.attr("dimnames");
  CharacterVector colNames = dimNames[1];
  CharacterVector rowNames = dimNames[0];
  int sizeHyperparam = hyperparam.ncol();
  CharacterVector sortedColNames(sizeHyperparam), sortedRowNames(sizeHyperparam);
  for(int i = 0; i < sizeHyperparam; i++)
    sortedColNames(i) = colNames(i), sortedRowNames(i) = rowNames(i);
  sortedColNames.sort();
  sortedRowNames.sort();
  
  for(int i = 0; i < sizeHyperparam; i++){
    if(i > 0 && (sortedColNames(i) == sortedColNames(i-1) || sortedRowNames(i) == sortedRowNames(i-1)))  
      stop("The states must all be unique");
    else if(sortedColNames(i) != sortedRowNames(i))
      stop("The set of row names must be the same as the set of column names");
    bool found = false;
    for(int j = 0; j < sizeMatr; j++)
      if(elements(j) == sortedColNames(i))
        found = true;
    // hyperparam may contain states not in stringchar
    if(!found)  elements.push_back(sortedColNames(i));
  }
  
  // check for the case where hyperparam has missing data
  for(int i = 0; i < sizeMatr; i++){
    bool found = false;
    for(int j = 0; j < sizeHyperparam; j++)
      if(sortedColNames(j) == elements(i))
        found = true;
    if(!found)
      stop("Hyperparameters for all state transitions must be provided");
  }   
      
  elements = elements.sort();
  sizeMatr = elements.size();
  
  for(int i = 0; i < sizeMatr; i++)
    for(int j = 0; j < sizeMatr; j++)
      if(hyperparam(i, j) < 1.)
        stop("The hyperparameter elements must all be greater than or equal to 1");
  
  // permute the elements of hyperparam such that the row, column names are sorted
  hyperparam = sortByDimNames(hyperparam);
  
  NumericMatrix freqMatr(sizeMatr), newFreqMatr(sizeMatr);

  double predictiveDist = 0.; // log of the predictive probability

  // populate frequeny matrix for old data; this is used for inference 
  int posFrom = 0, posTo = 0;
  for(int i = 0; i < stringchar.size() - 1; i ++) {
    for (int j = 0; j < sizeMatr; j ++) {
      if(stringchar[i] == elements[j]) posFrom = j;
      if(stringchar[i + 1] == elements[j]) posTo = j;
    }
    freqMatr(posFrom,posTo)++;
  }
  
  // frequency matrix for new data
  for(int i = 0; i < newData.size() - 1; i ++) {
    for (int j = 0; j < sizeMatr; j ++) {
      if(newData[i] == elements[j]) posFrom = j;
      if(newData[i + 1] == elements[j]) posTo = j;
    }
    newFreqMatr(posFrom,posTo)++;
  }
 
  for (int i = 0; i < sizeMatr; i++) {
    double rowSum = 0, newRowSum = 0, paramRowSum = 0;
    for (int j = 0; j < sizeMatr; j++){ 
      rowSum += freqMatr(i, j), newRowSum += newFreqMatr(i, j), paramRowSum += hyperparam(i, j);
      predictiveDist += lgamma(freqMatr(i, j) + newFreqMatr(i, j) + hyperparam(i, j)) -
                        lgamma(freqMatr(i, j) + hyperparam(i, j));
    }
    predictiveDist += lgamma(rowSum + paramRowSum) - lgamma(rowSum + newRowSum + paramRowSum);
  }

  return predictiveDist;
}


// [[Rcpp::export]]
NumericVector priorDistribution(NumericMatrix transMatr, NumericMatrix hyperparam = NumericMatrix()){
  // begin validity checks for the transition matrix
  if(transMatr.nrow() != transMatr.ncol())
    stop("Transition matrix dimensions are inconsistent");
    
  int sizeMatr = transMatr.nrow();
  for(int i = 0; i < sizeMatr; i++){
    double rowSum = 0., eps = 1e-10;
    for(int j = 0; j < sizeMatr; j++)
      if(transMatr(i, j) < 0. || transMatr(i, j) > 1.)
        stop("The entries in the transition matrix must each belong to the interval [0, 1]");
      else
        rowSum += transMatr(i, j);
    if(rowSum <= 1. - eps || rowSum >= 1. + eps)
      stop("The rows of the transition matrix must each sum to 1");
  }
  
  List dimNames = transMatr.attr("dimnames");
  if(dimNames.size() == 0)
    stop("Provide dimnames for the transition matrix");
  CharacterVector colNames = dimNames[1];
  CharacterVector rowNames = dimNames[0];
  CharacterVector sortedColNames(sizeMatr), sortedRowNames(sizeMatr);
  for(int i = 0; i < sizeMatr; i++)
    sortedColNames(i) = colNames(i), sortedRowNames(i) = rowNames(i);
  sortedColNames.sort();
  sortedRowNames.sort();
  
  for(int i = 0; i < sizeMatr; i++) 
    if(i > 0 && (sortedColNames(i) == sortedColNames(i-1) || sortedRowNames(i) == sortedRowNames(i-1)))  
      stop("The states must all be unique");
    else if(sortedColNames(i) != sortedRowNames(i))
      stop("The set of row names must be the same as the set of column names");
  
  // if no hyperparam argument provided, use default value of 1 for all 
  if(hyperparam.nrow() == 1 && hyperparam.ncol() == 1){
    NumericMatrix temp(sizeMatr, sizeMatr);
    temp.attr("dimnames") = List::create(sortedColNames, sortedColNames);
    for(int i = 0; i < sizeMatr; i++)
      for(int j = 0; j < sizeMatr; j++)
        temp(i, j) = 1;
    hyperparam = temp;
  }
  
  // validity check for hyperparam
  if(hyperparam.nrow() != hyperparam.ncol())
    stop("Dimensions of the hyperparameter matrix are inconsistent");
    
  if(hyperparam.nrow() != sizeMatr)
    stop("Hyperparameter and the transition matrices differ in dimensions");
    
  List _dimNames = hyperparam.attr("dimnames");
  if(_dimNames.size() == 0)
    stop("Provide dimnames for the hyperparameter matrix");
  CharacterVector _colNames = _dimNames[1];
  CharacterVector _rowNames = _dimNames[0];
  int sizeHyperparam = hyperparam.ncol();
  CharacterVector _sortedColNames(sizeHyperparam), _sortedRowNames(sizeHyperparam);
  for(int i = 0; i < sizeHyperparam; i++)
    _sortedColNames(i) = colNames(i), _sortedRowNames(i) = rowNames(i);
  _sortedColNames.sort();
  _sortedRowNames.sort();
  
  for(int i = 0; i < sizeHyperparam; i++)
    if(sortedColNames(i) != _sortedColNames(i) || sortedRowNames(i) != _sortedRowNames(i))
      stop("Hyperparameter and the transition matrices states differ");
  
  for(int i = 0; i < sizeMatr; i++)
    for(int j = 0; j < sizeMatr; j++)
      if(hyperparam(i, j) < 1.)
        stop("The hyperparameter elements must all be greater than or equal to 1");
 
  transMatr = sortByDimNames(transMatr);
  hyperparam = sortByDimNames(hyperparam);
  
  NumericVector logProbVec;
  for(int i = 0; i < sizeMatr; i++){
    double logProb_i = 0., hyperparamRowSum = 0;
    for(int j = 0; j < sizeMatr; j++){
      hyperparamRowSum += hyperparam(i, j);
      logProb_i += (hyperparam(i, j) - 1.) * log(transMatr(i, j)) - lgamma(hyperparam(i, j));
    }
    logProb_i += lgamma(hyperparamRowSum);
    logProbVec.push_back(logProb_i);
  }
  logProbVec.attr("names") = sortedColNames;

  return logProbVec;
}
