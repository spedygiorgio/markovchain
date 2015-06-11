List _mcFitMap(CharacterVector stringchar, bool byrow, double confidencelevel, NumericMatrix hyperparam) {
  // get mapEstMatr and freqMatr 
  CharacterVector elements = stringchar;
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
  std::sort(sortedColNames.begin(), sortedColNames.end());
  std::sort(sortedRowNames.begin(), sortedRowNames.end());
  
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
    
  // permute the elements of hyperparam such that the row, column names are sorted
  hyperparam = sortByDimNames(hyperparam);
  
  NumericMatrix mapEstMatr(sizeMatr), expMatr(sizeMatr);
  NumericMatrix freqMatr(sizeMatr);
  mapEstMatr.attr("dimnames") = List::create(elements, elements); 
  expMatr.attr("dimnames") = List::create(elements, elements); 

  NumericMatrix lowerEndpointMatr = NumericMatrix(mapEstMatr.nrow(), mapEstMatr.ncol());
  NumericMatrix upperEndpointMatr = NumericMatrix(mapEstMatr.nrow(), mapEstMatr.ncol());
  NumericMatrix varianceMatr = NumericMatrix(mapEstMatr.nrow(), mapEstMatr.ncol());

  // populate frequeny matrix for old data; this is used for inference 
  int posFrom = 0, posTo = 0;
  for(int i = 0; i < stringchar.size() - 1; i ++) {
    for (int j = 0; j < sizeMatr; j ++) {
      if(stringchar[i] == elements[j]) posFrom = j;
      if(stringchar[i + 1] == elements[j]) posTo = j;
    }
    freqMatr(posFrom,posTo)++;
  }
 
  // sanitize and to row probs
  for (int i = 0; i < sizeMatr; i++) {
    double rowSum = 0, paramRowSum = 0;
    for (int j = 0; j < sizeMatr; j++)
      rowSum += freqMatr(i, j), paramRowSum += hyperparam(i, j);
    
    // toRowProbs
    for (int j = 0; j < sizeMatr; j++) {
      // confidence intervals and bounds
      double p = freqMatr(i, j) + hyperparam(i, j), q = rowSum + paramRowSum - freqMatr(i, j) - hyperparam(i, j);
      
      // expected value of the transition parameters
      expMatr(i, j) = p / (p + q);
      
      if(p + q == sizeMatr)
              mapEstMatr(i, j) = 1 / sizeMatr;
      else
              // maximum a posteriori estimate
              mapEstMatr(i, j) = (p - 1) / (p + q - sizeMatr);
              
      double beta = lbeta(p, q);
      double cdf = betain(double(mapEstMatr(i, j)), p, q, beta);

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

  if(byrow==false) mapEstMatr = _transpose(mapEstMatr);

  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = mapEstMatr;
  outMc.slot("name") = "Bayesian Fit";  
  
  return List::create(_["estimate"] = outMc
    , _["confidenceInterval"] = List::create(_["confidenceLevel"]=confidencelevel, 
              _["lowerEndpointMatrix"]=lowerEndpointMatr, 
              _["upperEndpointMatrix"]=upperEndpointMatr),
              _["expectedValue"]=expMatr,
              _["varianceMatrix"]=varianceMatr
  );
}
