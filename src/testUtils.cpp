// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_set>
#include <string>
#include <iostream>
using namespace Rcpp;
using namespace arma;
using namespace std;

// This method receives the output of communicatingClasses(object) and object@states
// and checks that in fact the communicating classes are a partition of states
// Is a method agnostic on whether the Markov Chain was given by rows or columns
// [[Rcpp::export(.testthatIsPartitionRcpp)]]
bool isPartition(List commClasses, 
                 CharacterVector states) {
  int n = states.size();
  unordered_set<string> used;
  unordered_set<string> originalStates;
  int numClassStates = 0;
  bool partition = true;
  
  for (auto state : states)
    originalStates.insert((string) state);
  
  // Check that the union of the classes is
  // states and they do not overlap
  for (int i = 0; i < commClasses.size() && partition; ++i) {
    CharacterVector currentClass = commClasses(i);
    numClassStates += currentClass.size();

    for (int j = 0; j < currentClass.size() && partition; ++j) {
      string state = (string) currentClass(j);
      partition = used.count(state) == 0 && originalStates.count(state) > 0;
      used.insert(state);
    }
  }

  return partition && numClassStates == n;
}

// This is simply a method that checks the following recurrence,
// naming p = probs, f = hitting, it checks:
//
// f(i, j) = p(i, j) + ∑_{k ≠ j} p(i, k) f(k, j)
//
// where p are the transitionMatrix probs and hitting are the
// hitting probabilities for the Markov Chain associated to
// probs. byrow indicates whether probs is an stochastic matrix
// by rows or by columns. tolerance is the tolerance we want to
// check against (values must sufficiently close to 1)
// [[Rcpp::export(.testthatAreHittingRcpp)]]
bool areHittingProbabilities(NumericMatrix probs, 
                             NumericMatrix hitting,
                             bool byrow,
                             double tolerance) {
  if (!byrow) {
    probs = transpose(probs);
    hitting = transpose(hitting);
  }
  
  int numStates = probs.nrow();
  bool holds = true;
  double result;
  
  for (int i = 0; i < numStates && holds; ++i) {
    for (int j = 0; j < numStates && holds; ++j) {
      result = 0;
      
      for (int k = 0; k < numStates; ++k)
        if (k != j)
          result -= probs(i, k) * hitting(k, j);
        
      result += hitting(i, j) - probs(i, j);
      holds = abs(result) < tolerance;
    }
  }
  
  return holds;
}
