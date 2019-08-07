#include <Rcpp.h>
#include <ctime>

using namespace Rcpp;

#include <math.h>

List markovchainFit(SEXP data, String method = "mle", bool byrow = true,
                    int nboot = 10, double laplacian = 0, String name = "",
                    bool parallel = false, double confidencelevel = 0.95, bool confint = true,
                    NumericMatrix hyperparam = NumericMatrix(), bool sanitize = false,
                    CharacterVector possibleStates = CharacterVector()); 

//' @name ctmcFit
//' @title Function to fit a CTMC
//' @description This function fits the underlying CTMC give the state
//'   transition data and the transition times using the maximum likelihood
//'   method (MLE)
//' @usage ctmcFit(data, byrow = TRUE, name = "", confidencelevel = 0.95)
//' @param data It is a list of two elements. The first element is a character
//'   vector denoting the states. The second is a numeric vector denoting the
//'   corresponding transition times.
//' @param byRow Determines if the output transition probabilities of the
//'   underlying embedded DTMC are by row.
//' @param name Optional name for the CTMC.
//' @param confidencelevel Confidence level for the confidence interval
//'   construnction.
//' @return It returns a list containing the CTMC object and the confidence intervals.
//' 
//' @details  Note that in data, there must exist an element wise corresponding
//'   between the two elements of the list and that data[[2]][1] is always 0.
//' @references Continuous Time Markov Chains (vignette), Sai Bhargav Yalamanchi, Giorgio Alfredo Spedicato 2015
//' @author Sai Bhargav Yalamanchi
//' @seealso \code{\link{rctmc}}
//' 
//' @examples
//' data <- list(c("a", "b", "c", "a", "b", "a", "c", "b", "c"), c(0, 0.8, 2.1, 2.4, 4, 5, 5.9, 8.2, 9))
//' ctmcFit(data)
//' 
//' @export
//' 
// [[Rcpp::export]]
List ctmcFit(List data, bool byrow=true, String name="", double confidencelevel = 0.95) {
  
  CharacterVector stateData(as<CharacterVector>(data[0]).size());
  
  for (int i = 0; i < as<CharacterVector>(data[0]).size(); i++)
    stateData[i] = as<CharacterVector>(data[0])[i];

  NumericVector transData = data[1];
  CharacterVector sortedStates = unique(as<CharacterVector>(data[0])).sort();
  NumericVector stateCount(sortedStates.size());
  NumericVector stateSojournTime(sortedStates.size());
  
  List dtmcData = markovchainFit(stateData, "mle", byrow, 10, 0, name, false, confidencelevel);
  
  for (int i = 0; i < stateData.size() - 1; i++){
    int idx = std::find(sortedStates.begin(),
                        sortedStates.end(),
                        stateData[i]) - sortedStates.begin();
    stateCount[idx]++;
    stateSojournTime[idx] += transData[i+1] - transData[i];
  }
  
  S4 dtmcEst = dtmcData["estimate"];
  NumericMatrix gen = dtmcEst.slot("transitionMatrix");
  
  for (int i = 0; i < gen.nrow(); i++){
    for (int j = 0; j < gen.ncol(); j++){
      if (stateCount[i] > 0)
        gen(i, j) *= stateCount[i] / stateSojournTime[i];
    }
    if (stateCount[i] > 0)
      gen(i, i) = - stateCount[i] / stateSojournTime[i];
    else  
      gen(i, i) = -1;
  }
  
  double zscore = stats::qnorm_0(confidencelevel, 1.0, 0.0);
  NumericVector lowerConfVecLambda(sortedStates.size()), upperConfVecLambda(sortedStates.size());
  
  for (int i = 0; i < sortedStates.size(); i++){

    if (stateCount[i] > 0){
      auto factor = stateCount[i] / stateSojournTime[i] * (1 - zscore / sqrt(stateCount[i]));
      lowerConfVecLambda(i) = std::max(0., factor);
      upperConfVecLambda(i) = std::min(1., factor);
    } else {
      lowerConfVecLambda(i) = 1;
      upperConfVecLambda(i) = 1;
    }
  }
  
  S4 outCtmc("ctmc");
  outCtmc.slot("states") = sortedStates;
  outCtmc.slot("generator") = gen;
  outCtmc.slot("name") = name;
  
  return List::create(_["estimate"] = outCtmc,
                      _["errors"] = List::create(
                          _["dtmcConfidenceInterval"] = List::create(
                              _["confidenceLevel"] = dtmcData["confidenceLevel"],
                              _["lowerEndpointMatrix"] = dtmcData["lowerEndpointMatrix"],
                              _["upperEndpointMatrix"] = dtmcData["upperEndpointMatrix"]),
                          _["lambdaConfidenceInterval"] = List::create(
                              _["lowerEndpointVector"] = lowerConfVecLambda,
                              _["upperEndpointVector"] = upperConfVecLambda))
                      );
}
