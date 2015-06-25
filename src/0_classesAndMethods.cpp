// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export(.isProbRcpp)]]
bool isProb(double prob)
{
	if (prob<0 || prob >1) return false;
	return true;
}

// [[Rcpp::export]]
NumericMatrix generatorToTransitionMatrix(NumericMatrix gen){
  NumericMatrix transMatr(gen.nrow());
  transMatr.attr("dimnames") = gen.attr("dimnames");
  
  for(int i = 0; i < gen.nrow(); i++){
    for(int j = 0; j < gen.ncol(); j++){
      if(i != j)
        transMatr(i, j) = -gen(i, j) / gen(i, i);
    }
  }
  
  return transMatr;
}

// [[Rcpp::export(.isGenRcpp)]]
bool isGen(NumericMatrix gen){
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

/*// [[Rcpp::export]]
bool isValidTransitionMatrix(NumericMatrix transMatr){
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
  if(dimNames.size() == 2){
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
  }
 
  return true;
}

// [[Rcpp::export]]
bool isValidHyperparameterMatrix(NumericMatrix hyperparam, CharacterVector states = CharacterVector()){
  int sizeMatr = 0;
  if(states.size() > 0){
    states = unique(states).sort();
    sizeMatr = states.size();
  }
  
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
  }
  
  // check for the case where hyperparam has missing data
  for(int i = 0; i < sizeMatr; i++){
    bool found = false;
    for(int j = 0; j < sizeHyperparam; j++)
      if(sortedColNames(j) == states(i))
        found = true;
    if(!found)
      stop("Hyperparameters for all state transitions must be provided");
  }   
  
  for(int i = 0; i < sizeHyperparam; i++)
    for(int j = 0; j < sizeHyperparam; j++)
      if(hyperparam(i, j) < 1.)
        stop("The hyperparameter states must all be greater than or equal to 1");
        
  return true;
}*/