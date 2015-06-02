List _mcFitMap(CharacterVector stringchar, bool byrow, double confidencelevel, NumericMatrix hyperparam, CharacterVector newData) {
  // get initialMatr and freqMatr 
  CharacterVector elements = stringchar;
  for(int i = 0; i < newData.size(); i++)
    elements.push_back(newData[i]);
  
  elements = unique(elements).sort();
  int sizeMatr = elements.size();
  
  // if no hyperparam argument provided, use default value of 1 for all 
  if(hyperparam.nrow() == 1 && hyperparam.ncol() == 1){
    NumericMatrix temp(sizeMatr, sizeMatr);
    for(int i = 0; i < sizeMatr; i++)
      for(int j = 0; j < sizeMatr; j++)
        temp(i, j) = 1;
    hyperparam = temp;
  }
  
  // validity check
  if(hyperparam.nrow() != sizeMatr || hyperparam.ncol() != sizeMatr) 
    stop("Dimensions of the hyperparameter matrix are inconsistent");
  
  NumericMatrix initialMatr(sizeMatr);
  NumericMatrix freqMatr(sizeMatr), newFreqMatr(sizeMatr);
  initialMatr.attr("dimnames") = List::create(elements, elements); 

  NumericMatrix lowerEndpointMatr = NumericMatrix(initialMatr.nrow(), initialMatr.ncol());
  NumericMatrix upperEndpointMatr = NumericMatrix(initialMatr.nrow(), initialMatr.ncol());
  NumericMatrix varianceMatr = NumericMatrix(initialMatr.nrow(), initialMatr.ncol());
  
  double predictiveDist = 0.; // log of the predictive probability

  // populate frequeny matrix for old data; this is used for inference 
  int posFrom, posTo;
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
 
  // sanitize and to row probs
  for (int i = 0; i < sizeMatr; i++) {
    double rowSum = 0, newRowSum = 0, paramRowSum = 0;
    for (int j = 0; j < sizeMatr; j++){ 
      rowSum += freqMatr(i, j), newRowSum += newFreqMatr(i, j), paramRowSum += hyperparam(i, j);
      predictiveDist += lgamma(freqMatr(i, j) + newFreqMatr(i, j) + hyperparam(i, j)) -
                        lgamma(freqMatr(i, j) + hyperparam(i, j));
    }
    
    predictiveDist += lgamma(rowSum + paramRowSum) - lgamma(rowSum + newRowSum + paramRowSum);
    // toRowProbs
    for (int j = 0; j < sizeMatr; j++) {
      if(rowSum == 0)
              initialMatr(i, j) = 1 / sizeMatr;
      else
              // maximum a posteriori estimate
              initialMatr(i, j) = (freqMatr(i, j) + hyperparam(i, j)) / (rowSum + paramRowSum);

      // confidence intervals and bounds
      double p = freqMatr(i, j) + hyperparam(i, j), q = rowSum + paramRowSum - freqMatr(i, j) - hyperparam(i, j);
      double beta = lbeta(p, q);
      double cdf = betain(double(initialMatr(i, j)), p, q, beta);

      if(cdf + confidencelevel / 2 > 1.){
        upperEndpointMatr(i, j) = 1.;
        lowerEndpointMatr(i, j) = xinbta(p, q, beta, 1 - confidencelevel);
      }
      else if(cdf - confidencelevel / 2 < 0.){
        lowerEndpointMatr(i, j) = 0.;
        upperEndpointMatr(i, j) = xinbta(p, q, beta, confidencelevel);
      }
      else{
        lowerEndpointMatr(i, j) = xinbta(p, q, beta, cdf - confidencelevel / 2);
        upperEndpointMatr(i, j) = xinbta(p, q, beta, cdf + confidencelevel / 2);
      }

      varianceMatr(i, j) = p * q / (p + q) / (p + q) / (1 + p + q);
    }
  }

  if(byrow==false) initialMatr = _transpose(initialMatr);

  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = initialMatr;
  outMc.slot("name") = "Bayesian Fit";  
  
  return List::create(_["estimate"] = outMc
    , _["confidenceInterval"] = List::create(_["confidenceLevel"]=confidencelevel, 
              _["lowerEndpointMatrix"]=lowerEndpointMatr, 
              _["upperEndpointMatrix"]=upperEndpointMatr),
              _["varianceMatrix"]=varianceMatr,
              _["predictiveProbability"]=exp(predictiveDist)
  );
}
