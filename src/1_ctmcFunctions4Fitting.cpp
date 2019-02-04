#include <Rcpp.h>
#include <ctime>

using namespace Rcpp;

#include <math.h>

List markovchainFit(SEXP data, String method="mle", bool byrow=true, int nboot=10, double laplacian=0
            , String name="", bool parallel=false, double confidencelevel=0.95, bool confint = true
            , NumericMatrix hyperparam = NumericMatrix(), bool sanitize = false, CharacterVector possibleStates = CharacterVector()); 

// [[Rcpp::export]]
List ctmcFit(List data, bool byrow=true, String name="", double confidencelevel = 0.95)
{
  CharacterVector stateData(as<CharacterVector>(data[0]).size());
  for(int i = 0; i < as<CharacterVector>(data[0]).size(); i++)
    stateData[i] = as<CharacterVector>(data[0])[i];
  NumericVector transData = data[1];
  CharacterVector sortedStates = unique(as<CharacterVector>(data[0])).sort();
  NumericVector stateCount(sortedStates.size());
  NumericVector stateSojournTime(sortedStates.size());
  
  List dtmcData = markovchainFit(stateData, "mle", byrow, 10, 0, name, false, confidencelevel);
  
  for(int i = 0; i < stateData.size() - 1; i++){
    int idx = std::find(sortedStates.begin(), sortedStates.end(), stateData[i]) - sortedStates.begin();
    stateCount[idx]++;
    stateSojournTime[idx] += transData[i+1] - transData[i];
  }
  
  S4 dtmcEst = dtmcData["estimate"];
  NumericMatrix gen = dtmcEst.slot("transitionMatrix");
  
  for(int i = 0; i < gen.nrow(); i++){
    for(int j = 0; j < gen.ncol(); j++){
      if(stateCount[i] > 0)
        gen(i, j) *= stateCount[i] / stateSojournTime[i];
    }
    if(stateCount[i] > 0)
      gen(i, i) = - stateCount[i] / stateSojournTime[i];
    else  
      gen(i, i) = -1;
  }
  
  double zscore = stats::qnorm_0(confidencelevel, 1.0, 0.0);
  NumericVector lowerConfVecLambda(sortedStates.size()), upperConfVecLambda(sortedStates.size());
  
  for(int i = 0; i < sortedStates.size(); i++){
    if(stateCount[i] > 0){
      lowerConfVecLambda(i) = std::max(0., stateCount[i] / stateSojournTime[i] * (1 - zscore / sqrt(stateCount[i])));
      upperConfVecLambda(i) = std::min(1., stateCount[i] / stateSojournTime[i] * (1 + zscore / sqrt(stateCount[i])));
    }
    else{
      lowerConfVecLambda(i) = 1;
      upperConfVecLambda(i) = 1;
    }
  }
  
  S4 outCtmc("ctmc");
  outCtmc.slot("states") = sortedStates;
  outCtmc.slot("generator") = gen;
  outCtmc.slot("name") = name;
  
  return List::create(_["estimate"] = outCtmc,
                      _["errors"] = List::create(_["dtmcConfidenceInterval"] = List::create(
                                                 _["confidenceLevel"] = dtmcData["confidenceLevel"],
                                                 _["lowerEndpointMatrix"] = dtmcData["lowerEndpointMatrix"],
                                                 _["upperEndpointMatrix"] = dtmcData["upperEndpointMatrix"]),
                      _["lambdaConfidenceInterval"] = List::create(_["lowerEndpointVector"] = lowerConfVecLambda,
                      _["upperEndpointVector"] = upperConfVecLambda)));
}