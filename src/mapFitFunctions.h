/*
 * Function to remove NA values from a vector
 */

CharacterVector clean_nas(CharacterVector elements_na){
  int count = 0;
  // first cycle: count the number of non-NA elements
  for(int i = 0; i < elements_na.size(); i++) {
    if(elements_na[i] != "NA") count++;
  }
  
  // create a new vector with the correct size
  CharacterVector elements(count);
  int idx = 0;
  // Second cycle: populate the vector
  for(int i = 0; i < elements_na.size(); i++) {
    if(elements_na[i] != "NA") {
      elements[idx++] = elements_na[i];
    }
  }
  return elements;
}
  

List _mcFitMap(SEXP data, bool byrow, double confidencelevel, NumericMatrix hyperparam = NumericMatrix(), 
               bool sanitize = false, CharacterVector possibleStates = CharacterVector()) {
  
  if(TYPEOF(data) != VECSXP)  {
    data = List::create(as<CharacterVector>(data));
  }
  
  List seqs = as<List>(data);
  CharacterVector elements;
  
  for(int i = 0;i < (int)seqs.size();i++) {
    CharacterVector tseq = unique(as<CharacterVector>(seqs[i]));
    for(int j = 0;j < (int)tseq.size();j++) {
      if(tseq[j] != "NA") {
        elements.push_back(tseq[j]);
      }
    }
  }

  elements = unique(union_(elements, possibleStates)).sort();
  // number of unique states
  int sizeMatr = elements.size();
  
  // if no hyperparam argument provided, use default value of 1 for all.
  // The R-level default is `matrix()`, which is a 1x1 matrix containing
  // NA -- checking the NA is essential here, not just the 1x1 shape:
  // otherwise a genuinely user-supplied 1x1 hyperparam matrix (e.g. one
  // that mistakenly omits states actually present in the data, which
  // happens to leave it 1x1) would be silently discarded and replaced by
  // the default instead of being validated and rejected below.
  if(hyperparam.nrow() == 1 && hyperparam.ncol() == 1 && NumericVector::is_na(hyperparam(0, 0))) {
    // matrix with all entries 1
    NumericMatrix temp(sizeMatr, sizeMatr);
    temp.attr("dimnames") = List::create(elements, elements);
    for(int i = 0; i < sizeMatr; i++)
      for(int j = 0; j < sizeMatr; j++)
        temp(i, j) = 1;
    hyperparam = temp;
  }
  
  //-----------beginning of validity checking of hyperparam matrix----------------------
  
  // validity check for hyperparam matrix
  if(hyperparam.nrow() != hyperparam.ncol()) {
    stop("Dimensions of the hyperparameter matrix are inconsistent");
  }
    
  if(hyperparam.nrow() < sizeMatr) {
    stop("Hyperparameters for all state transitions must be provided");
  }
    
  // extract rows and columns name out of hyperparam matrix   
  List dimNames = hyperparam.attr("dimnames");
  CharacterVector colNames = dimNames[1];
  CharacterVector rowNames = dimNames[0];
  
  // size of hyperparam matrix
  int sizeHyperparam = hyperparam.ncol();
  
  // sorted order of hyperparam rows and columns name
  CharacterVector sortedColNames(sizeHyperparam), sortedRowNames(sizeHyperparam);
  for(int i = 0; i < sizeHyperparam; i++) {
    sortedColNames(i) = colNames(i), sortedRowNames(i) = rowNames(i);
  }
  sortedColNames.sort();
  sortedRowNames.sort();
  
  // validity of hyperparam matrix
  for(int i = 0; i < sizeHyperparam; i++){
    
    // columns names must be different
    // rows names must be different
    if(i > 0 && (sortedColNames(i) == sortedColNames(i-1) || sortedRowNames(i) == sortedRowNames(i-1))) {
      stop("The states must all be unique");
    } 
    
    // same states should be present in rows and columns   
    else if(sortedColNames(i) != sortedRowNames(i)) {
      stop("The set of row names must be the same as the set of column names");
    }
      
    // chech whether any state in column names exists which is not in the given sequence
    bool found = false;
    for(int j = 0; j < sizeMatr; j++) {
      if(elements(j) == sortedColNames(i)) {
        found = true;
        break;
      }
    }
      
    // hyperparam may contain states not in stringchar
    if(!found) {
      elements.push_back(sortedColNames(i)); 
    }
  }
  
  // check for the case where hyperparam has missing data
  for(int i = 0; i < sizeMatr; i++){
    bool found = false;
    for(int j = 0; j < sizeHyperparam; j++) {
      if(sortedColNames(j) == elements(i)) {
        found = true;
        break;
      }
    }
      
    if(!found)
      stop("Hyperparameters for all state transitions must be provided");
  }   
      
  
  elements = elements.sort();
  sizeMatr = elements.size();
  
  for(int i = 0; i < sizeMatr; i++)
    for(int j = 0; j < sizeMatr; j++)
      if(hyperparam(i, j) < 1.)
        stop("The hyperparameter elements must all be greater than or equal to 1");
    
  //-----------end of validity checking of hyperparam matrix----------------------
  
  // permute the elements of hyperparam such that the row, column names are sorted
  hyperparam = sortByDimNames(hyperparam);
  
  // helper matrices which will help in the calculation of transition matrix
  // other matrices will be returned as a result
  NumericMatrix mapEstMatr(sizeMatr), expMatr(sizeMatr);
  NumericMatrix freqMatr(sizeMatr);
  mapEstMatr.attr("dimnames") = List::create(elements, elements); 
  expMatr.attr("dimnames") = List::create(elements, elements); 

  // matrices to be returned
  NumericMatrix lowerEndpointMatr = NumericMatrix(mapEstMatr.nrow(), mapEstMatr.ncol());
  NumericMatrix upperEndpointMatr = NumericMatrix(mapEstMatr.nrow(), mapEstMatr.ncol());
  NumericMatrix stdError = NumericMatrix(mapEstMatr.nrow(), mapEstMatr.ncol());

  // populate frequeny matrix for old data; this is used for inference
  for(int k = 0;k < seqs.size();k++) {
    CharacterVector stringchar  = as<CharacterVector>(seqs[k]);
    int posFrom = 0, posTo = 0;
    for(R_xlen_t i = 0; i < stringchar.size() - 1; i ++) {
      if(stringchar[i] != "NA" && stringchar[i+1] != "NA"){
        for (int j = 0; j < sizeMatr; j ++) {
          if(stringchar[i] == elements[j]) posFrom = j;
          if(stringchar[i + 1] == elements[j]) posTo = j;
        }
        freqMatr(posFrom,posTo)++;
      }
    }  
  }
 
  // sanitize and to row probs
  for (int i = 0; i < sizeMatr; i++) {
    
    // rowsum of frequency matrix and hyperparam matrix
    double rowSum = 0, paramRowSum = 0;
    for (int j = 0; j < sizeMatr; j++) {
      rowSum += freqMatr(i, j), paramRowSum += hyperparam(i, j);
    }
      
    // toRowProbs
    for (int j = 0; j < sizeMatr; j++) {
      
      // confidence intervals and bounds
      double p = freqMatr(i, j) + hyperparam(i, j), q = rowSum + paramRowSum - freqMatr(i, j) - hyperparam(i, j);
      
      // expected value of the transition parameters
      expMatr(i, j) = p / (p + q);
      
      if(p + q == sizeMatr) {
        mapEstMatr(i, j) = (sanitize ? 1.0 / sizeMatr : 0);
      }
      else {
        // maximum a posteriori estimate
        mapEstMatr(i, j) = (p - 1) / (p + q - sizeMatr);
      }

      // Equal-tailed credible interval: quantiles of the exact marginal
      // Beta(p, q) posterior at (1-confidencelevel)/2 and
      // 1-(1-confidencelevel)/2. (The previous construction centered the
      // interval, in CDF space, at the posterior MODE rather than the
      // median; that has exact nominal coverage too, since it is still of
      // the form [F^-1(u - cl/2), F^-1(u + cl/2)] for a fixed u, but it is
      // non-standard, needs special-casing near the boundary, and -- 
      // verified by simulation -- tends to produce wider intervals with
      // *worse* empirical coverage than the equal-tailed interval below.)
      double beta = lbeta(p, q);
      double halfAlpha = (1.0 - confidencelevel) / 2.0;
      lowerEndpointMatr(i, j) = xinbta(p, q, beta, halfAlpha);
      upperEndpointMatr(i, j) = xinbta(p, q, beta, 1.0 - halfAlpha);

      stdError(i, j) = sqrt(p * q / (p + q) / (p + q) / (1 + p + q));
    }
  }

  lowerEndpointMatr.attr("dimnames") = upperEndpointMatr.attr("dimnames")
    = stdError.attr("dimnames") = List::create(elements, elements);

  // transpose the matrix if columwise result is required
  if(byrow == false) {
    mapEstMatr = transposeMatrix(mapEstMatr); 
  }

  // markovchain object
  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = mapEstMatr;
  outMc.slot("name") = "Bayesian Fit";  
  
 // message("\n\'estimate\' is the MAP set of parameters where as \'expectedValue\' \nis the expectation 
 // of the parameters with respect to the posterior.\nThe confidence intervals are given for \'estimate\'.");
  
  return List::create(_["estimate"] = outMc,
                      _["expectedValue"] = expMatr,
                      _["standardError"] = stdError,
                      _["confidenceInterval"] = List::create(_["confidenceLevel"] = confidencelevel, 
                                                             _["lowerEndpointMatrix"] = lowerEndpointMatr, 
                                                             _["upperEndpointMatrix"] = upperEndpointMatr));
}
