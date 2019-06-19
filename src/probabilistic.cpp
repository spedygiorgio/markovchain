// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <math.h>
#include <unordered_set>
#include <string>
#include <algorithm>
#include <stack>

using namespace Rcpp;
using namespace std;

template <typename T>
T sortByDimNames(const T m);

// check if two vectors are intersected
bool _intersected(CharacterVector x, CharacterVector y) {
  if (x.size() < y.size())
    return _intersected(y, x);
  else {
    unordered_set<string> values;
    bool intersect = false;
    
    for (auto value : x)
      values.insert(as<string>(value));
   
    for (auto it = y.begin(); it != y.end() && !intersect; ++it)
      intersect = values.count(as<string>(*it)) > 0;
  
    return intersect;
  }
}

// [[Rcpp::export(.commClassesKernelRcpp)]]
SEXP commClassesKernel(NumericMatrix P) {
  unsigned int numStates = P.ncol();
  CharacterVector stateNames = rownames(P);
  int numReachable;
  int classSize;
  
  // The entry (i,j) of this matrix is true iff we can reach j from i
  vector<vector<bool>> communicates(numStates, vector<bool>(numStates, false));
  vector<list<int>> adjacencies(numStates);
  
  // We fill the adjacencies matrix for the graph
  // A state j is in the adjacency of i iff P(i, j) > 0
  for (int i = 0; i < numStates; ++i)
    for (int j = 0; j < numStates; ++j)
      if (P(i, j) > 0)
        adjacencies[i].push_back(j);


  // Backtrack from all the states to find which
  // states communicate with a given on
  // O(n³) where n is the number of states
  for (int i = 0; i < numStates; ++i) {
    stack<int> notVisited;
    notVisited.push(i);
    
    while (!notVisited.empty()) {
      int j = notVisited.top();
      notVisited.pop();
      communicates[i][j] = true;
      
      for (int k: adjacencies[j])
        if (!communicates[i][k])
          notVisited.push(k);
    }
  }
  
  LogicalMatrix classes(numStates, numStates);
  classes.attr("dimnames") = List::create(stateNames, stateNames);
  // v populated with FALSEs
  LogicalVector closed(numStates);
  closed.names() = stateNames;
  
  for (int i = 0; i < numStates; ++i) {
    numReachable = 0;
    classSize = 0;
    
    /* We mark i and j as the same communicating class iff we can reach the
       state j from i and the state i from j
       We count the size of the communicating class of i (i is fixed here),
       and if it matches the number of states that can be reached from i,
       then the class is closed
    */
    for (int j = 0; j < numStates; ++j) {
      classes(i, j) = communicates[i][j] && communicates[j][i];
      
      if (classes(i,j))
        classSize += 1;

      // Number of states reachable from i
      if (communicates[i][j])
        numReachable += 1;
    }
    
    if (classSize == numReachable)
      closed(i) = true;
  }
  
  return List::create(_["classes"] = classes, _["closed"] = closed);
}

// [[Rcpp::export(.communicatingClassesRcpp)]]
List communicatingClasses(S4 object) {
  //returns the underlying communicating classes
  
  NumericMatrix matr = object.slot("transitionMatrix");
  List temp = commClassesKernel(matr);
  LogicalMatrix adjMatr = temp["classes"];
  int len = adjMatr.nrow();
  List classesList;
  CharacterVector rnames = rownames(adjMatr);
  
  for (int i = 0; i < len; i ++) {
    bool isNull = false;
    LogicalVector row2Check = adjMatr(i, _);
    CharacterVector proposedCommClass;
    
    for (int j = 0; j < row2Check.size(); j++) {
      if (row2Check[j] == true) {
        String rname = rnames[j];
        proposedCommClass.push_back(rname);
      }
    }
    
    if (i > 0) {
      for (int j = 0; j < classesList.size(); j ++) {
        bool check = false;        
        CharacterVector cv = classesList[j];
        std::set<std::string> s1, s2;
        
        for (int k = 0; k < cv.size(); k ++) {
          s1.insert(as<std::string>(cv[k]));
          if (proposedCommClass.size() > k) {
            s2.insert(as<std::string>(proposedCommClass[k]));
          }
        }
        
        check = std::equal(s1.begin(), s1.end(), s2.begin());
        
        if (check) {
          isNull = true;
          break;
        }
      }
    }
    
    if (!isNull) 
      classesList.push_back(proposedCommClass);
  }
  return classesList;
}

// returns the recurrent classes
// [[Rcpp::export(.recurrentClassesRcpp)]]
List recurrentClasses(S4 object) {
  NumericMatrix matr = object.slot("transitionMatrix");
  List temp = commClassesKernel(matr);
  List communicatingClassList = communicatingClasses(object);
  List v = temp["closed"];
  CharacterVector ns = v.names();
  CharacterVector transientStates;
  
  for (int i = 0; i < v.size(); i++) {
    if (bool(v[i]) == false)
      transientStates.push_back(ns[i]);
  }
  
  List recurrentClassesList;
  
  for (int i = 0; i < communicatingClassList.size(); i ++) {
    CharacterVector class2Test = communicatingClassList[i];
    
    if (!_intersected(class2Test,transientStates)) 
      recurrentClassesList.push_back(class2Test);
  }
  
  return recurrentClassesList;
}

// matrix power function
arma::mat _pow(arma::mat A, int n) {
  arma::mat R = arma::eye(A.n_rows, A.n_rows);
  
  for (int i = 0; i < n; i ++) 
    R = A*R;
  
  return R;
}

//communicating states
// [[Rcpp::export(.commStatesFinderRcpp)]]
NumericMatrix commStatesFinder(NumericMatrix matr) {
  //Reachability matrix
  int dimMatr = matr.nrow();
  arma::mat X(matr.begin(), dimMatr, dimMatr, false);
  arma::mat temp = arma::eye(dimMatr, dimMatr) + arma::sign(X);
  temp = _pow(temp, dimMatr - 1);
  NumericMatrix R = wrap(arma::sign(temp));
  R.attr("dimnames") = matr.attr("dimnames");
  
  return R;
}

// summary of markovian object
// [[Rcpp::export(.summaryKernelRcpp)]]
List summaryKernel(S4 object) {
  NumericMatrix matr = object.slot("transitionMatrix");
  List temp = commClassesKernel(matr);
  List communicatingClassList = communicatingClasses(object);
  // List communicatingClassList = communicatingClasses(temp["C"]);
  List v = temp["closed"];
  CharacterVector ns = v.names();
  CharacterVector transientStates;
  
  for (int i = 0; i < v.size(); i++) {
    if (bool(v[i]) == false)
      transientStates.push_back(ns[i]);
  }
  List closedClasses, transientClasses, recurrentClassesList;

  for (int i = 0; i < communicatingClassList.size(); i ++) {
    CharacterVector class2Test = communicatingClassList[i];
    
    if (_intersected(class2Test,transientStates)) 
        transientClasses.push_back(class2Test); 
      else {
        bool isClosed = true;
        
        for (int j = 0; j < class2Test.size(); j ++) {
          for (int from = 0; from < ns.size(); from ++) {
            bool inclass = false;
            
            if (ns[from] != class2Test[j]) continue;
            
            for (int to = 0; to < matr.cols(); to ++) {
              for (int k = 0; k < class2Test.size(); k ++) {
                if (class2Test[k] == ns[to]) inclass = true;
              }
              if (from == to || inclass) continue;
              if (matr(from, to) != 0) {
                isClosed = false;
                break;
              }
            }
          }
        }
        if (isClosed)
          closedClasses.push_back(class2Test);
        // recurrentClassesList.push_back(class2Test);
      }
  }
  recurrentClassesList = recurrentClasses(object);
  List summaryMc = List::create(_["closedClasses"] = closedClasses,
                                _["recurrentClasses"] = recurrentClassesList,
                                _["transientClasses"] = transientClasses);
  
  return(summaryMc);
}

//here the kernel function to compute the first passage
// [[Rcpp::export(.firstpassageKernelRcpp)]]
NumericMatrix firstpassageKernel(NumericMatrix P, int i, int n) {
  arma::mat G = as<arma::mat>(P);
  arma::mat Pa = G;
  arma::mat H(n, P.ncol()); 
  
  //here Thoralf suggestion
  //initializing the first row
  for (unsigned int j = 0; j < G.n_cols; j++)
    H(0, j) = G(i-1, j);
  
  arma::mat E = 1 - arma::eye(P.ncol(), P.ncol());

  for (int m = 1; m < n; m++) {
    G = Pa * (G%E);
    
    for (unsigned int j = 0; j < G.n_cols; j ++) 
      H(m, j) = G(i-1, j);
  }
  
  NumericMatrix R = wrap(H);
  
  return R;
}



// [[Rcpp::export(.firstPassageMultipleRCpp)]]
NumericVector firstPassageMultipleRCpp(NumericMatrix P,int i, NumericVector setno, int n) {
  arma::mat G = as<arma::mat>(P);
  arma::mat Pa = G;
  arma::vec H = arma::zeros(n); //here Thoralf suggestion
  unsigned int size = setno.size();
  //initializing the first row
  for (unsigned int k = 0; k < size; k++) {
    H[0] += G(i-1, setno[k]-1);
  }
  
  arma::mat E = 1 - arma::eye(P.ncol(), P.ncol());
  
  for (int m = 1; m < n; m++) {
    G = Pa * (G%E);
    
    for (unsigned int k = 0; k < size; k++) {
      H[m] += G(i-1, setno[k]-1);
    }
  }
  
  NumericVector R = wrap(H);
  
  return R;
}

// [[Rcpp::export(.expectedRewardsRCpp)]]
NumericVector expectedRewardsRCpp(NumericMatrix matrix, int n, NumericVector rewards) {
  // initialises output vector
  NumericVector out;
  
  // gets no of states
  int no_of_states = matrix.ncol();
  
  // initialises armadillo matrices and vectors
  arma::vec temp = arma::zeros(no_of_states);
  arma::mat matr = as<arma::mat>(matrix);
  arma::vec v = arma::zeros(no_of_states);
  
  // initialses the vector for the base case of dynamic programming expression
  for (int i=0;i<no_of_states;i++) {
    temp[i] = rewards[i];
    v[i] = rewards[i];
  }
  
  // v(n, u) = r + [P]v(n−1, u);
  for (int i=0;i<n;i++) {
    temp = v + matr*temp;
  }
  
  // gets output in form of NumericVector
  out = wrap(temp);
  
  return out;
}

// [[Rcpp::export(.expectedRewardsBeforeHittingARCpp)]]
double expectedRewardsBeforeHittingARCpp(NumericMatrix matrix,int s0,
                               NumericVector rewards, int n ) {
  float result = 0.0;
  int size = rewards.size();
  arma::mat matr = as<arma::mat>(matrix);
  arma::mat temp = as<arma::mat>(matrix);
  arma::vec r = as<arma::vec>(rewards);
  arma::mat I = arma::zeros(1,size);
  
  I(0,s0-1) = 1;
  
  for (int j = 0; j < n; j++) {
    arma::mat res = I*(temp*r);
    result = result + res(0,0);
    temp = temp*matr;
  }
  
  return result;
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

// function to get the period of a DTMC

//' @rdname absorbingStates
//' 
//' @export
//' 
// [[Rcpp::export(period)]]
int period(S4 object) {
  Function isIrreducible("is.irreducible");
  List res = isIrreducible(object);
  
  if (!res[0]) {
    warning("The matrix is not irreducible");
    return 0;
  } else {
    NumericMatrix P = object.slot("transitionMatrix");
    int n = P.ncol();
    std::vector<double> r, T(1), w;
    int d = 0, m = T.size(), i = 0, j = 0;
    
    if (n > 0) {
      arma::vec v(n);
      v[0] = 1;
      
      while (m>0 && d!=1) {
        i = T[0];
        T.erase(T.begin());
        w.push_back(i);
        j = 0;
        
        while (j < n) {
          if (P(i,j) > 0) {
            r.insert(r.end(), w.begin(), w.end());
            r.insert(r.end(), T.begin(), T.end());
            double k = 0;
            
            for (std::vector<double>::iterator it = r.begin(); it != r.end(); it ++) 
              if (*it == j) k ++;
            
            if (k > 0) {
               int b = v[i] + 1 - v[j];
               d = gcd(d, b);
            } else {
              T.push_back(j);
              v[j] = v[i] + 1;
            }
          }
          j++;
        }
        m = T.size();
      }
    }
    
    // v = v - floor(v/d)*d;
    return d;
  }
}

// [[Rcpp::export]]
double predictiveDistribution(CharacterVector stringchar, CharacterVector newData, NumericMatrix hyperparam = NumericMatrix()) {
  // construct list of states
  CharacterVector elements = stringchar;
  
  for (int i = 0; i < newData.size(); i++)
    elements.push_back(newData[i]);
  
  elements = unique(elements).sort();
  int sizeMatr = elements.size();
  
  // if no hyperparam argument provided, use default value of 1 for all 
  if (hyperparam.nrow() == 1 && hyperparam.ncol() == 1) {
    NumericMatrix temp(sizeMatr, sizeMatr);
    temp.attr("dimnames") = List::create(elements, elements);
    
    for (int i = 0; i < sizeMatr; i++)
      for (int j = 0; j < sizeMatr; j++)
        temp(i, j) = 1;
    
    hyperparam = temp;
  }
  
  // validity check
  if (hyperparam.nrow() != hyperparam.ncol())
    stop("Dimensions of the hyperparameter matrix are inconsistent");
    
  if (hyperparam.nrow() < sizeMatr)
    stop("Hyperparameters for all state transitions must be provided");
    
  List dimNames = hyperparam.attr("dimnames");
  CharacterVector colNames = dimNames[1];
  CharacterVector rowNames = dimNames[0];
  int sizeHyperparam = hyperparam.ncol();
  CharacterVector sortedColNames(sizeHyperparam), sortedRowNames(sizeHyperparam);
  
  for (int i = 0; i < sizeHyperparam; i++)
    sortedColNames(i) = colNames(i), sortedRowNames(i) = rowNames(i);

  sortedColNames.sort();
  sortedRowNames.sort();
  
  for (int i = 0; i < sizeHyperparam; i++) {
    if (i > 0 && (sortedColNames(i) == sortedColNames(i-1) || sortedRowNames(i) == sortedRowNames(i-1)))  
      stop("The states must all be unique");
    else if (sortedColNames(i) != sortedRowNames(i))
      stop("The set of row names must be the same as the set of column names");
    
    bool found = false;
    
    for (int j = 0; j < sizeMatr; j++)
      if (elements(j) == sortedColNames(i))
        found = true;
    // hyperparam may contain states not in stringchar
    if (!found)  elements.push_back(sortedColNames(i));
  }
  
  // check for the case where hyperparam has missing data
  for (int i = 0; i < sizeMatr; i++) {
    bool found = false;
    
    for (int j = 0; j < sizeHyperparam; j++)
      if (sortedColNames(j) == elements(i))
        found = true;
    
    if (!found)
      stop("Hyperparameters for all state transitions must be provided");
  }   
      
  elements = elements.sort();
  sizeMatr = elements.size();
  
  for (int i = 0; i < sizeMatr; i++)
    for (int j = 0; j < sizeMatr; j++)
      if (hyperparam(i, j) < 1.)
        stop("The hyperparameter elements must all be greater than or equal to 1");
  
  // permute the elements of hyperparam such that the row, column names are sorted
  hyperparam = sortByDimNames(hyperparam);
  
  NumericMatrix freqMatr(sizeMatr), newFreqMatr(sizeMatr);

  double predictiveDist = 0.; // log of the predictive probability

  // populate frequeny matrix for old data; this is used for inference 
  int posFrom = 0, posTo = 0;
  
  for (int i = 0; i < stringchar.size() - 1; i ++) {
    for (int j = 0; j < sizeMatr; j ++) {
      if (stringchar[i] == elements[j]) posFrom = j;
      if (stringchar[i + 1] == elements[j]) posTo = j;
    }
    freqMatr(posFrom,posTo)++;
  }
  
  // frequency matrix for new data
  for (int i = 0; i < newData.size() - 1; i ++) {
    for (int j = 0; j < sizeMatr; j ++) {
      if (newData[i] == elements[j]) posFrom = j;
      if (newData[i + 1] == elements[j]) posTo = j;
    }
    newFreqMatr(posFrom,posTo)++;
  }
 
  for (int i = 0; i < sizeMatr; i++) {
    double rowSum = 0, newRowSum = 0, paramRowSum = 0;
    
    for (int j = 0; j < sizeMatr; j++) { 
      rowSum += freqMatr(i, j), newRowSum += newFreqMatr(i, j), paramRowSum += hyperparam(i, j);
      predictiveDist += lgamma(freqMatr(i, j) + newFreqMatr(i, j) + hyperparam(i, j)) -
                        lgamma(freqMatr(i, j) + hyperparam(i, j));
    }
    predictiveDist += lgamma(rowSum + paramRowSum) - lgamma(rowSum + newRowSum + paramRowSum);
  }

  return predictiveDist;
}


// [[Rcpp::export]]
NumericVector priorDistribution(NumericMatrix transMatr, NumericMatrix hyperparam = NumericMatrix()) {
  // begin validity checks for the transition matrix
  if (transMatr.nrow() != transMatr.ncol())
    stop("Transition matrix dimensions are inconsistent");
    
  int sizeMatr = transMatr.nrow();
  
  for (int i = 0; i < sizeMatr; i++) {
    double rowSum = 0., eps = 1e-10;
    
    for (int j = 0; j < sizeMatr; j++)
      if (transMatr(i, j) < 0. || transMatr(i, j) > 1.)
        stop("The entries in the transition matrix must each belong to the interval [0, 1]");
      else
        rowSum += transMatr(i, j);
    
    if (rowSum <= 1. - eps || rowSum >= 1. + eps)
      stop("The rows of the transition matrix must each sum to 1");
  }
  
  List dimNames = transMatr.attr("dimnames");
  
  if (dimNames.size() == 0)
    stop("Provide dimnames for the transition matrix");
  
  CharacterVector colNames = dimNames[1];
  CharacterVector rowNames = dimNames[0];
  CharacterVector sortedColNames(sizeMatr), sortedRowNames(sizeMatr);
  
  for (int i = 0; i < sizeMatr; i++)
    sortedColNames(i) = colNames(i), sortedRowNames(i) = rowNames(i);
  
  sortedColNames.sort();
  sortedRowNames.sort();
  
  for (int i = 0; i < sizeMatr; i++) 
    if (i > 0 && (sortedColNames(i) == sortedColNames(i-1) || sortedRowNames(i) == sortedRowNames(i-1)))  
      stop("The states must all be unique");
    else if (sortedColNames(i) != sortedRowNames(i))
      stop("The set of row names must be the same as the set of column names");
  
  // if no hyperparam argument provided, use default value of 1 for all 
  if (hyperparam.nrow() == 1 && hyperparam.ncol() == 1) {
    NumericMatrix temp(sizeMatr, sizeMatr);
    temp.attr("dimnames") = List::create(sortedColNames, sortedColNames);
  
    for (int i = 0; i < sizeMatr; i++)
      for (int j = 0; j < sizeMatr; j++)
        temp(i, j) = 1;
  
    hyperparam = temp;
  }
  
  // validity check for hyperparam
  if (hyperparam.nrow() != hyperparam.ncol())
    stop("Dimensions of the hyperparameter matrix are inconsistent");
    
  if (hyperparam.nrow() != sizeMatr)
    stop("Hyperparameter and the transition matrices differ in dimensions");
    
  List _dimNames = hyperparam.attr("dimnames");

  if (_dimNames.size() == 0)
    stop("Provide dimnames for the hyperparameter matrix");
  
  CharacterVector _colNames = _dimNames[1];
  CharacterVector _rowNames = _dimNames[0];
  int sizeHyperparam = hyperparam.ncol();
  CharacterVector _sortedColNames(sizeHyperparam), _sortedRowNames(sizeHyperparam);
  
  for (int i = 0; i < sizeHyperparam; i++)
    _sortedColNames(i) = colNames(i), _sortedRowNames(i) = rowNames(i);
  
  _sortedColNames.sort();
  _sortedRowNames.sort();
  
  for (int i = 0; i < sizeHyperparam; i++)
    if (sortedColNames(i) != _sortedColNames(i) || sortedRowNames(i) != _sortedRowNames(i))
      stop("Hyperparameter and the transition matrices states differ");
  
  for (int i = 0; i < sizeMatr; i++)
    for (int j = 0; j < sizeMatr; j++)
      if (hyperparam(i, j) < 1.)
        stop("The hyperparameter elements must all be greater than or equal to 1");
 
  transMatr = sortByDimNames(transMatr);
  hyperparam = sortByDimNames(hyperparam);
  NumericVector logProbVec;
  
  for (int i = 0; i < sizeMatr; i++) {
    double logProb_i = 0., hyperparamRowSum = 0;
  
    for (int j = 0; j < sizeMatr; j++) {
      hyperparamRowSum += hyperparam(i, j);
      logProb_i += (hyperparam(i, j) - 1.) * log(transMatr(i, j)) - lgamma(hyperparam(i, j));
    }
    
    logProb_i += lgamma(hyperparamRowSum);
    logProbVec.push_back(logProb_i);
  }
  
  logProbVec.attr("names") = sortedColNames;

  return logProbVec;
}

// [[Rcpp::export(.hittingProbabilitiesRcpp)]]
NumericMatrix hittingProbabilities(NumericMatrix transitionMatrix) {
  int numStates = transitionMatrix.nrow();
  arma::mat transitionProbs = as<arma::mat>(transitionMatrix);
  arma::mat result(numStates, numStates);
  // Compute closed communicating classes
  List commClasses = commClassesKernel(transitionMatrix);
  List closedClass = commClasses["closed"];
  LogicalMatrix communicating = commClasses["classes"];

  
  for (int j = 0; j < numStates; ++j) {
    arma::mat to_invert = as<arma::mat>(transitionMatrix);
    arma::vec right_part = -transitionProbs.col(j);
    
    for (int i = 0; i < numStates; ++i) {
      to_invert(i, j) = 0;
      to_invert(i, i) -= 1;
    }

    for (int i = 0; i < numStates; ++i) {
      if (closedClass(i)) {
        for (int k = 0; k < numStates; ++k)
          if (k != i)
            to_invert(i, k) = 0;
          else
            to_invert(i, i) = 1;
          
        if (communicating(i, j))
          right_part(i) = 1;
        else
          right_part(i) = 0;
      }
    }
    
    arma::mat inverse = arma::inv(to_invert);
    result.col(j) = inverse * right_part;
  }
  
  return wrap(result);
}
