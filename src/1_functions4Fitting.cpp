// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <Rcpp.h>

using namespace Rcpp;

#include "mathHelperFunctions.h"
#include "mapFitFunctionsSAI.h"
#include <math.h>

NumericMatrix _toRowProbs(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericMatrix out(nrow);

  for (int i = 0; i < nrow; i++) {
    double rowSum = 0;
    for (int j = 0; j < ncol; j++) 
      rowSum += x(i, j);
    for (int j = 0; j < ncol; j++) 
      out(i, j) = x(i, j)/rowSum;
  }
  out.attr("dimnames") = List::create(rownames(x), colnames(x)); 
  return out;
}

// [[Rcpp::export]]
NumericMatrix createSequenceMatrix(CharacterVector stringchar, bool toRowProbs=false, bool sanitize=true) {
  CharacterVector elements = unique(stringchar).sort();
  int sizeMatr = elements.size();
  
  NumericMatrix freqMatrix(sizeMatr);
  freqMatrix.attr("dimnames") = List::create(elements, elements); 
  CharacterVector rnames = rownames(freqMatrix);

  int posFrom = 0, posTo = 0;
  for(int i = 0; i < stringchar.size() - 1; i ++) {
    for (int j = 0; j < rnames.size(); j ++) {
      if(stringchar[i] == rnames[j]) posFrom = j;
      if(stringchar[i + 1] == rnames[j]) posTo = j;
    }
  	freqMatrix(posFrom,posTo)++;
  }
 
  //sanitizing if any row in the matrix sums to zero by posing the corresponding diagonal equal to 1/dim
  if(sanitize==true)
  {
    for (int i = 0; i < sizeMatr; i++) {
      double rowSum = 0;
      for (int j = 0; j < sizeMatr; j++) 
      rowSum += freqMatrix(i, j);
      if(rowSum == 0)
      for (int j = 0; j < sizeMatr; j++) 
      freqMatrix(i, j) = 1/sizeMatr;
    }
  }
  if(toRowProbs==true)
  freqMatrix = _toRowProbs(freqMatrix);
  
  return (freqMatrix);
}

double _loglikelihood(CharacterVector seq, NumericMatrix transMatr) {
  double out = 0;
  CharacterVector rnames = rownames(transMatr);
  int from = 0, to = 0; 
  for(int i = 0; i < seq.size() - 1; i ++) {
    for(int r = 0; r < rnames.size(); r ++) {
      if(rnames[r] == seq[i]) from = r; 
      if(rnames[r] == seq[i + 1]) to = r; 
    }    
    out += log(transMatr(from, to));
  }
  return out;
}

List _mcFitMle(CharacterVector stringchar, bool byrow, double confidencelevel) {
  // get initialMatr and freqMatr 
  CharacterVector elements = unique(stringchar).sort();
  int sizeMatr = elements.size();
  
  NumericMatrix initialMatr(sizeMatr);
  NumericMatrix freqMatr(sizeMatr);
  initialMatr.attr("dimnames") = List::create(elements, elements); 

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
  	double rowSum = 0;
  	for (int j = 0; j < sizeMatr; j++) 
  		rowSum += freqMatr(i, j);
  	// toRowProbs
  	for (int j = 0; j < sizeMatr; j++) {
  	  if(rowSum == 0)
  	    initialMatr(i, j) = 1/sizeMatr;
  	  else
  	    initialMatr(i, j) = freqMatr(i, j)/rowSum;
  	}
  }

  if(byrow==false) initialMatr = _transpose(initialMatr);

  // confidence interval
  double zscore = stats::qnorm_0(confidencelevel, 1.0, 0.0);
  int nrows = initialMatr.nrow();
  int ncols = initialMatr.ncol();
  NumericMatrix lowerEndpointMatr(nrows, ncols);
  NumericMatrix upperEndpointMatr(nrows, ncols);
  NumericMatrix standardError(nrows, ncols);

  double marginOfError, lowerEndpoint, upperEndpoint;
  for(int i = 0; i < nrows; i ++) {
    for(int j = 0; j < ncols; j ++) {
      standardError(i, j) = initialMatr(i, j) / sqrt(freqMatr(i, j));
      marginOfError = zscore * standardError(i, j);
      lowerEndpoint = initialMatr(i, j) - marginOfError;
      upperEndpoint = initialMatr(i, j) + marginOfError;
      lowerEndpointMatr(i,j) = (lowerEndpoint > 1.0) ? 1.0 : ((0.0 > lowerEndpoint) ? 0.0 : lowerEndpoint);
      upperEndpointMatr(i,j) = (upperEndpoint > 1.0) ? 1.0 : ((0.0 > upperEndpoint) ? 0.0 : upperEndpoint);
    }
  }
  standardError.attr("dimnames") = upperEndpointMatr.attr("dimnames") 
          = lowerEndpointMatr.attr("dimnames") = List::create(elements, elements); 

  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = initialMatr;
  outMc.slot("name") = "MLE Fit";  
  
  return List::create(_["estimate"] = outMc
          , _["standardError"] = standardError
		      , _["confidenceInterval"] = List::create(_["confidenceLevel"]=confidencelevel, 
		        _["lowerEndpointMatrix"]=lowerEndpointMatr, 
		        _["upperEndpointMatrix"]=upperEndpointMatr)							
	);
}

List _mcFitLaplacianSmooth(CharacterVector stringchar, bool byrow, double laplacian=0.01) {
  NumericMatrix origNum = createSequenceMatrix(stringchar, false);
  int nRows = origNum.nrow(), nCols = origNum.ncol();
  for(int i = 0; i < nRows; i ++) {
	double rowSum = 0;
	for(int j = 0; j < nCols; j ++) {
    		origNum(i,j) += laplacian;
    		rowSum += origNum(i,j);
  	}
  	//get a transition matrix and a DTMC
	for(int j = 0; j < nCols; j ++) 
    		origNum(i,j) /= rowSum;
  }
  
  if(byrow==false) origNum = _transpose(origNum);
 
  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = origNum;
  outMc.slot("name") = "Laplacian Smooth Fit";  

  return List::create(_["estimate"] = outMc);
}

List _bootstrapCharacterSequences(CharacterVector stringchar, int n, int size=-1) {
  if(size == -1) size = stringchar.size();
  NumericMatrix contingencyMatrix = createSequenceMatrix(stringchar);
  List samples, res;
  CharacterVector itemset = rownames(contingencyMatrix);
  int itemsetsize = itemset.size();

  Function sample("sample");
  for(int i = 0; i < n; i ++) {
    CharacterVector charseq, resvec;	
    int rnd = (int)(runif(1)(0) * itemsetsize);
    String ch = itemset[rnd];
    charseq.push_back(ch);
    for(int j = 1; j < size; j ++) {
      NumericVector probsVector;
      for(int k = 0; k < itemsetsize; k ++) {
        if((std::string)itemset[k] == (std::string) ch) {
          probsVector = contingencyMatrix(k, _);	
          break;
        }
      }
      res = sample(itemset, 1, true, probsVector);
      resvec = res[0];
      ch = resvec[0];
      charseq.push_back(ch);
    }
    samples.push_back(charseq);
  }

  return samples;
}

List _fromBoot2Estimate(List listMatr) {
  int sampleSize = listMatr.size();
  NumericMatrix firstMat = listMatr[0];
  int matrDim = firstMat.nrow();
  NumericMatrix matrMean(matrDim), matrSd(matrDim);

  for(int i = 0; i < matrDim; i ++) { 
  	for(int j = 0; j < matrDim; j ++) { 
  	  NumericVector probsEstimated;
  	  for(int k = 0; k < sampleSize; k ++) {
  	    NumericMatrix mat = listMatr[k];
  	    probsEstimated.push_back(mat(i,j));
  	  }
  	  matrMean(i,j) = mean(probsEstimated);
  	  matrSd(i,j) = sd(probsEstimated);
  	}
  }
  matrMean.attr("dimnames") = List::create(rownames(firstMat), colnames(firstMat)); 
  matrSd.attr("dimnames") = matrMean.attr("dimnames");
  return List::create(_["estMu"]=matrMean, _["estSigma"]=matrSd);
}

// worker for parallel loop
struct ForLoopWorker : public RcppParallel::Worker
{
   const List input;
   List output;

   ForLoopWorker(const List input, List output)
      : input(input), output(output) {}

   void operator()(std::size_t begin, std::size_t end) {
	output[begin] = createSequenceMatrix(input[begin], true, true);
   }
};

List _mcFitBootStrap(CharacterVector data, int nboot, bool byrow, bool parallel, double confidencelevel) {
  List theList = _bootstrapCharacterSequences(data, nboot);
  int n = theList.size();
  List pmsBootStrapped(n);

  if(!parallel) { 
    for(int i = 0; i < n; i++) 
      pmsBootStrapped[i] = createSequenceMatrix(theList[i], true, true);
  } else {
  	ForLoopWorker forloop(theList, pmsBootStrapped);
  	parallelFor(0, n, forloop);
  }
  List estimateList = _fromBoot2Estimate(pmsBootStrapped);
  NumericMatrix transMatr = _toRowProbs(estimateList["estMu"]);

  S4 estimate("markovchain");
  estimate.slot("transitionMatrix") = transMatr;
  estimate.slot("byrow") = byrow;
  estimate.slot("name") = "BootStrap Estimate";  

  // confidence interval
  double zscore = stats::qnorm_0(confidencelevel, 1.0, 0.0);
  int nrows = transMatr.nrow();
  int ncols = transMatr.ncol();
  NumericMatrix lowerEndpointMatr(nrows, ncols), upperEndpointMatr(nrows, ncols);
  NumericMatrix sigma = estimateList["estSigma"], standardError(nrows, ncols);
  
  double marginOfError, lowerEndpoint, upperEndpoint;
  for(int i = 0; i < nrows; i ++) {
    for(int j = 0; j < ncols; j ++) {
      standardError(i, j) = sigma(i, j) / sqrt(n); 
      marginOfError = zscore * standardError(i, j);
      lowerEndpoint = transMatr(i, j) - marginOfError;
      upperEndpoint = transMatr(i, j) + marginOfError;
      lowerEndpointMatr(i,j) = (lowerEndpoint > 1.0) ? 1.0 : ((0.0 > lowerEndpoint) ? 0.0 : lowerEndpoint);
      upperEndpointMatr(i,j) = (upperEndpoint > 1.0) ? 1.0 : ((0.0 > upperEndpoint) ? 0.0 : upperEndpoint);
    }
  }
  standardError.attr("dimnames") = upperEndpointMatr.attr("dimnames") 
              = lowerEndpointMatr.attr("dimnames") = transMatr.attr("dimnames"); 
              
  return List::create(_["estimate"] = estimate
    , _["standardError"] = standardError
    , _["confidenceInterval"] = List::create(_["confidenceLevel"]=confidencelevel, 
      _["lowerEndpointMatrix"]=lowerEndpointMatr
      , _["upperEndpointMatrix"]=upperEndpointMatr)
    , _["bootStrapSamples"] = pmsBootStrapped
		);
}

S4 _matr2Mc(CharacterMatrix matrData, double laplacian=0) {
  int nRows = matrData.nrow(), nCols = matrData.ncol();

  std::set<std::string> uniqueVals;
  for(int i = 0; i < nRows; i++) 
  	for(int j = 0; j < nCols; j++) 
		uniqueVals.insert((std::string)matrData(i, j));	

  int usize = uniqueVals.size();
  NumericMatrix contingencyMatrix (usize);
  contingencyMatrix.attr("dimnames") = List::create(uniqueVals, uniqueVals); 
  
  std::set<std::string>::iterator it;
  int stateBegin = 0, stateEnd = 0;
  for(int i = 0; i < nRows; i ++) {
    for(int j = 1; j < nCols; j ++) {
      int k = 0;
      for(it=uniqueVals.begin(); it!=uniqueVals.end(); ++it, k++) {
        if(*it == (std::string)matrData(i,j-1)) stateBegin = k;
        if(*it == (std::string)matrData(i,j)) stateEnd = k;
      }
      contingencyMatrix(stateBegin,stateEnd)++;
    }
  }

  //add laplacian correction if needed
  for(int i = 0; i < usize; i ++) {
	double rowSum = 0;
	for(int j = 0; j < usize; j ++) {
    		contingencyMatrix(i,j) += laplacian;
    		rowSum += contingencyMatrix(i,j);
  	}
  	//get a transition matrix and a DTMC
	for(int j = 0; j < usize; j ++) 
    		contingencyMatrix(i,j) /= rowSum;
  }
  
  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = contingencyMatrix;

  return(outMc);
}

// [[Rcpp::export]]
List inferHyperparam(NumericMatrix transMatr = NumericMatrix(), NumericVector scale = NumericVector(), CharacterVector data = CharacterVector()){
  if(transMatr.nrow() * transMatr.ncol() == 1 && data.size() == 0)
    stop("Provide the prior transition matrix or the data set in order to infer the hyperparameters");
  
  List out;
  
  if(transMatr.nrow() * transMatr.ncol() != 1){
    if(scale.size() == 0)
      stop("Provide a non-zero scaling factor vector to infer integer hyperparameters");
      
    // begin validity checks for the transition matrix
    if(transMatr.nrow() != transMatr.ncol())
      stop("Transition matrix dimensions are inconsistent");
      
    int sizeMatr = transMatr.nrow();
    for(int i = 0; i < sizeMatr; i++){
      double rowSum = 0., eps = 1e-10;
      for(int j = 0; j < sizeMatr; j++)
        rowSum += transMatr(i, j);
      if(rowSum <= 1. - eps || rowSum >= 1. + eps)
        stop("The rows of the transition matrix must each sum to 1");
    }
    
    List dimNames = transMatr.attr("dimnames");
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
        
    // validity check for the scaling factor vector
    if(scale.size() != sizeMatr)
      stop("The dimensions of the scale vector must match the number of states in the chain");
      
    for(int i = 0; i < sizeMatr; i++)
      if(scale(i) == 0)
        stop("The scaling factors must be non-zero!");
    
    NumericMatrix hpScaled(sizeMatr);
    hpScaled.attr("dimnames") = List::create(rowNames, colNames);
    for(int i = 0; i < sizeMatr; i++)
      for(int j = 0; j < sizeMatr; j++)
        hpScaled(i, j) = scale(i) * transMatr(i, j);
        
    hpScaled = sortByDimNames(hpScaled);
    
    out = List::create(_["scaledInference"] = hpScaled);
  }
  
  else if(data.size() != 0){
    CharacterVector elements = data;
    for(int i = 0; i < data.size(); i++)
      elements.push_back(data[i]);
    
    elements = unique(elements).sort();
    int sizeMatr = elements.size();
    
    NumericMatrix hpData(sizeMatr);
    hpData.attr("dimnames") = List::create(elements, elements); 
    std::fill(hpData.begin(), hpData.end(), 1);
    
    int posFrom = 0, posTo = 0;
    for(int i = 0; i < data.size() - 1; i ++) {
      for (int j = 0; j < sizeMatr; j ++) {
        if(data[i] == elements[j]) posFrom = j;
        if(data[i + 1] == elements[j]) posTo = j;
      }
      hpData(posFrom,posTo)++;
    }
    
    out = List::create(_["dataInference"] = hpData);
  }
  
  return out;
}

// [[Rcpp::export]]
List markovchainFit(SEXP data, String method="mle", bool byrow=true, int nboot=10, double laplacian=0
            , String name="", bool parallel=false, double confidencelevel=0.95
            , NumericMatrix hyperparam = NumericMatrix()) 
{
  List out;
  if(Rf_inherits(data, "data.frame") || Rf_inherits(data, "matrix")) { 
    CharacterMatrix mat;
    //if data is a data.frame forced to matrix
  	if(Rf_inherits(data, "data.frame")) {
  	  DataFrame df(data);
  	  mat = CharacterMatrix(df.nrows(), df.size());
  	  for(int i = 0; i < df.size(); i++)
    	  mat(_,i) = CharacterVector(df[i]);
 	  } else {
		  mat = data;
	  }
  	//byrow assumes distinct observations (trajectiories) are per row
  	//otherwise transpose
  	if(!byrow) mat = _transpose(mat);
   	S4 outMc =_matr2Mc(mat,laplacian);
   	out = List::create(_["estimate"] = outMc);
  } else {
    if(method == "mle") out = _mcFitMle(data, byrow, confidencelevel);
    if(method == "bootstrap") out = _mcFitBootStrap(data, nboot, byrow, parallel, confidencelevel);
    if(method == "laplace") out = _mcFitLaplacianSmooth(data, byrow, laplacian);
    if(method == "map") out = _mcFitMap(data, byrow, confidencelevel, hyperparam);
  }

  if(name != "") {
    S4 estimate = out["estimate"];
    estimate.slot("name") = name;
    out["estimate"] = estimate;
  }
  
  S4 estimate = out["estimate"];
  NumericMatrix transMatr = estimate.slot("transitionMatrix");
  estimate.slot("states") = rownames(transMatr);
  out["estimate"] = estimate;
  if(!Rf_inherits(data, "data.frame") && !Rf_inherits(data, "matrix")) 
    out["logLikelihood"] = _loglikelihood(data, transMatr);
  return out;
}
