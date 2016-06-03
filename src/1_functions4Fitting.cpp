// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <ctime>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace RcppArmadillo;

#include "mathHelperFunctions.h"
#include "mapFitFunctionsSAI.h"
#include <math.h>


// [[Rcpp::export(.markovchainSequenceRcpp)]]
CharacterVector markovchainSequenceRcpp(int n, S4 markovchain, CharacterVector t0,
                                        bool include_t0 = false) {
  
  // character vector to store the result
  CharacterVector chain(n);
  
  // transition mastrix
  NumericMatrix transitionMatrix = markovchain.slot("transitionMatrix");
  
  // possible states
  CharacterVector states = markovchain.slot("states");
  
  // current state
  CharacterVector state = t0;
  
  NumericVector rowProbs(states.size());
  CharacterVector outstate;
  
  
  
  for(int i = 0;i < n;i++) {
    
    // extracting row probabilties for the given state from transition matrix
    int row_no = 0;
    for(int j = 0;j < states.size();j++) {
      if(states[j] == state[0]) {
        row_no = j;
        break;
      }
    }
    
    for(int j = 0; j < states.size(); j++) {
      rowProbs[j] = transitionMatrix(row_no, j);
    }
    
    // calculate next state
    outstate = sample(states, 1, false, rowProbs);
    chain[i] = outstate[0];  
    state = outstate;
    
  }
  
  if (include_t0)
    chain.push_front(t0[0]);
  
  return chain;
}

bool checkSequenceRcpp(List object) {
  bool out = true;
  int nob = object.size();
  
  // if there is only one markovchain object return true
  if (nob == 1)
    return(true);
  
  S4 ob0, ob1;
  CharacterVector statesNm1, statesN, intersection;
  
  for(int i = 1; i < nob;i++) {
    ob0 = S4(object[i-1]);
    ob1 = S4(object[i]);
    
    statesNm1 = ob0.slot("states"); 
    statesN = ob1.slot("states");
    
    intersection = intersect(statesNm1, statesN);
    if(not setequal(intersection, statesNm1)) {
      out = false;
      break;
    }
  }
  return(out);
}

// [[Rcpp::export(.markovchainListRcpp)]]
List markovchainListRcpp(int n, List object, bool include_t0 = false) {
  bool verify = checkSequenceRcpp(object);
  
  if (not verify) {
    warning("Warning: some states in the markovchain sequences are not contained in the following states!");
  }
  
  
  NumericVector iteration = NumericVector::create();
  CharacterVector values = CharacterVector::create();
  S4 ob(object[0]);
  
  CharacterVector sampledValues, newVals;
  IntegerVector outIter;
  
  for(int i = 0;i < n;i++) {
    sampledValues = markovchainSequenceRcpp(1, object[0], CharacterVector(ob.slot("states")), include_t0);
    outIter = rep(i+1, sampledValues.size());
    
    if(object.size() > 1) {
      for(int j = 1;j < object.size();j++) {
        newVals = markovchainSequenceRcpp(1, object[j], sampledValues, include_t0);
        outIter.push_back(i+1);
        sampledValues.push_back(newVals[0]);
      }
    }
    
    for(int k = 0;k < outIter.size();k++) {
      iteration.push_back(outIter[k]);
      values.push_back(sampledValues[k]);
    }
  }
  
  return(List::create(iteration, values));
}

// convert a frequency matrix to a transition probability matrix
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

// Create a frequency matrix
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
            freqMatrix(i, j) = 1;
      }
  }
  if(toRowProbs==true)
    freqMatrix = _toRowProbs(freqMatrix);
  
  return (freqMatrix);
}

// log-likelihood
double _loglikelihood(CharacterVector seq, NumericMatrix transMatr) {
  
  // to store the result
  double out = 0;
  
  // states names
  CharacterVector rnames = rownames(transMatr);
  
  // caculate out
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

// Fit DTMC using MLE
List _mcFitMle(CharacterVector stringchar, bool byrow, double confidencelevel) {
  
  // unique states
  CharacterVector elements = unique(stringchar).sort();
  // matrix size = nrows = ncols
  int sizeMatr = elements.size(); 
  
  // initial matrix = transition matrix
  NumericMatrix initialMatr(sizeMatr);
  // frequencey matrix 
  NumericMatrix freqMatr(sizeMatr);
  
  // set names of states as rows name and columns name
  initialMatr.attr("dimnames") = List::create(elements, elements); 

  // populate frequency matrix
  int posFrom = 0, posTo = 0; 
  for(int i = 0; i < stringchar.size() - 1; i++) {  // ### size ###
    for (int j = 0; j < sizeMatr; j++) {            
      if(stringchar[i] == elements[j]) posFrom = j;
      if(stringchar[i + 1] == elements[j]) posTo = j;
    }
    
    freqMatr(posFrom, posTo)++;
  }

  // take care of rows with all entries 0
  for (int i = 0; i < sizeMatr; i++) {  
  	double rowSum = 0;
  	for (int j = 0; j < sizeMatr; j++) { 
  	  rowSum += freqMatr(i, j);
  	}
  		
  	// calculate rows probability
  	for (int j = 0; j < sizeMatr; j++) {  
  	  if(rowSum == 0) {
  	    initialMatr(i, j) = 1.0/sizeMatr;
  	  }
  	  else {
  	    initialMatr(i, j) = freqMatr(i, j)/rowSum;
  	  }
  	}
  }

  // transpose the matrix if byrow is false
  if(byrow == false) {
    initialMatr = _transpose(initialMatr); 
  }

  // z score for given confidence interval
  double zscore = stats::qnorm_0(confidencelevel, 1.0, 0.0);
  
  // store dimension of matrix
  int nrows = initialMatr.nrow();
  int ncols = initialMatr.ncol();
  
  // matrix to store end results
  NumericMatrix lowerEndpointMatr(nrows, ncols);
  NumericMatrix upperEndpointMatr(nrows, ncols);
  NumericMatrix standardError(nrows, ncols);

  // populate above defined matrix 
  double marginOfError, lowerEndpoint, upperEndpoint;
  for(int i = 0; i < nrows; i++) {
    for(int j = 0; j < ncols; j++) {
      if(freqMatr(i, j) == 0) {
        
        // whether entire ith row is zero or not
        bool notrans = true;
        for(int k = 0; k < ncols; k++) {
          
          // if entire ith row is not zero then set notrans to false  
          if(freqMatr(i, k) != 0) {
            standardError(i, j) = lowerEndpointMatr(i, j) = upperEndpointMatr(i, j) = 0;
            notrans = false;
            break;
          }
          
        }
        
        // if entire ith row is zero 
        if(notrans) 
          standardError(i, j) = lowerEndpointMatr(i, j) = upperEndpointMatr(i, j) = 1;
      } 
      else {
        // standard error calculation
        standardError(i, j) = initialMatr(i, j) / sqrt(freqMatr(i, j));
        
        // marginal error calculation
        marginOfError = zscore * standardError(i, j);
        
        // lower and upper end point calculation
        lowerEndpoint = initialMatr(i, j) - marginOfError;
        upperEndpoint = initialMatr(i, j) + marginOfError;
        
        // taking care that upper and lower end point should be between 0(included) and 1(included)
        lowerEndpointMatr(i, j) = (lowerEndpoint > 1.0) ? 1.0 : ((0.0 > lowerEndpoint) ? 0.0 : lowerEndpoint);
        upperEndpointMatr(i, j) = (upperEndpoint > 1.0) ? 1.0 : ((0.0 > upperEndpoint) ? 0.0 : upperEndpoint);
      }
    }
  }
  
  // set the rows and columns name as states names
  standardError.attr("dimnames") = upperEndpointMatr.attr("dimnames") 
          = lowerEndpointMatr.attr("dimnames") = List::create(elements, elements); 

  // create markov chain object
  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = initialMatr;
  outMc.slot("name") = "MLE Fit";  
  
  // return a list of important results
  return List::create(_["estimate"] = outMc,
                      _["standardError"] = standardError,
		                  _["confidenceInterval"] = List::create(_["confidenceLevel"] = confidencelevel, 
		                                                         _["lowerEndpointMatrix"] = lowerEndpointMatr, 
		                                                         _["upperEndpointMatrix"] = upperEndpointMatr)							
	       );
}

// Fit DTMC using Laplacian smooth
List _mcFitLaplacianSmooth(CharacterVector stringchar, bool byrow, double laplacian = 0.01) {
  
  // create frequency matrix
  NumericMatrix origNum = createSequenceMatrix(stringchar, false);
  
  // store dimension of frequenct matrix
  int nRows = origNum.nrow(), nCols = origNum.ncol();
  
  // convert frequency matrix to transition matrix
  for(int i = 0; i < nRows; i ++) {
	  double rowSum = 0;
	  
	  // add laplacian correction to each entry
	  // also calculate row's sum
	  for(int j = 0; j < nCols; j ++) {
    		origNum(i,j) += laplacian;
    		rowSum += origNum(i,j);
  	}
	  
  	// get a transition matrix and a DTMC
	  for(int j = 0; j < nCols; j ++) {
	    origNum(i,j) /= rowSum; 
	  }
  }
  
  // transpose transition matrix = columnwise storage 
  if(byrow == false) {
    origNum = _transpose(origNum);  
  }
 
  // create markovchain object
  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = origNum;
  outMc.slot("name") = "Laplacian Smooth Fit";  

  return List::create(_["estimate"] = outMc);
}

// bootstrap a sequence to produce a list of sample sequences
List _bootstrapCharacterSequences(CharacterVector stringchar, int n, int size = -1) {
  
  // store length of sequence
  if(size == -1) {
    size = stringchar.size();  
  }
  
  // frequency matrix
  NumericMatrix contingencyMatrix = createSequenceMatrix(stringchar);
  
  // many samples from a given a sequence :: bootstrap
  // res list is helper list
  List samples, res;
  
  // state names
  CharacterVector itemset = rownames(contingencyMatrix);
  
  // number of distinct states
  int itemsetsize = itemset.size();

  // access R sample function
  Function sample("sample");
  
  
  for(int i = 0; i < n; i ++) {
    // charseq to store a fresh sequence
    CharacterVector charseq, resvec;	
    
    // select a random number between 0 ans itemsize
    int rnd = (int)(runif(1)(0) * itemsetsize);
    
    // random state
    String ch = itemset[rnd];
    
    // push start state to charseq
    charseq.push_back(ch);
    
    
    for(int j = 1; j < size; j ++) {
      // store row probability
      NumericVector probsVector;
      
      // populate probsVector
      for(int k = 0; k < itemsetsize; k ++) {
        if((std::string)itemset[k] == (std::string) ch) {
          probsVector = contingencyMatrix(k, _);	
          break;
        }
      }
      
      // select next state of sequence
      res = sample(itemset, 1, true, probsVector);
      resvec = res[0];
      
      // current state
      ch = resvec[0];
      charseq.push_back(ch);
    }
    
    // every add one sequence
    samples.push_back(charseq);
  }

  // return a list of n sequence of same length as of given sequence
  return samples;
}


// estimate from the list of bootstrapped matrices
List _fromBoot2Estimate(List listMatr) {
  
  // number of transition matrices 
  int sampleSize = listMatr.size();
  
  // first transition matrix
  NumericMatrix firstMat = listMatr[0];
  
  // dimension of matrix
  int matrDim = firstMat.nrow();
  
  // matrix to store mean and standard deviation
  NumericMatrix matrMean(matrDim), matrSd(matrDim);

  // populate mean and sd matrix
  for(int i = 0; i < matrDim; i ++) { 
  	for(int j = 0; j < matrDim; j ++) { 
  	  
  	  NumericVector probsEstimated;
  	  for(int k = 0; k < sampleSize; k ++) {
        NumericMatrix mat = listMatr[k];
  	    probsEstimated.push_back(mat(i,j));
  	  }
  	  
  	  matrMean(i, j) = mean(probsEstimated);
  	  matrSd(i, j) = Rcpp::sd(probsEstimated);
  	}
  }
  
  // set rows and columns names = states names
  matrMean.attr("dimnames") = List::create(rownames(firstMat), colnames(firstMat)); 
  matrSd.attr("dimnames") = matrMean.attr("dimnames");
  
  // return list of estimated mean transition matrix and estimated sd matrix 
  return List::create(_["estMu"] = matrMean, _["estSigma"] = matrSd);
}

// Fit DTMC using bootstrap method
List _mcFitBootStrap(CharacterVector data, int nboot, bool byrow, bool parallel, double confidencelevel) {
  
  // list of sequence generated using given sequence
  List theList = _bootstrapCharacterSequences(data, nboot);
  
  // number of new sequence
  int n = theList.size();
  
  // to store frequency matrix for every sequence
  List pmsBootStrapped(n);

  // populate pmsBootStrapped
  if(parallel)
    for(int i = 0; i < n; i++)
      pmsBootStrapped[i] = createSequenceMatrix(theList[i], true, true);
  
  else 
    for(int i = 0; i < n; i++) 
      pmsBootStrapped[i] = createSequenceMatrix(theList[i], true, true);
  
  
  List estimateList = _fromBoot2Estimate(pmsBootStrapped);
  
  // transition matrix
  NumericMatrix transMatr = _toRowProbs(estimateList["estMu"]);

  // markovchain object
  S4 estimate("markovchain");
  estimate.slot("transitionMatrix") = transMatr;
  estimate.slot("byrow") = byrow;
  estimate.slot("name") = "BootStrap Estimate";  

  // z score for given confidence interval
  double zscore = stats::qnorm_0(confidencelevel, 1.0, 0.0);
  
  // store dimension of matrix
  int nrows = transMatr.nrow();
  int ncols = transMatr.ncol();
  
  // matrix to store end results
  NumericMatrix lowerEndpointMatr(nrows, ncols), upperEndpointMatr(nrows, ncols);
  NumericMatrix sigma = estimateList["estSigma"], standardError(nrows, ncols);
  
  // populate above defined matrix 
  double marginOfError, lowerEndpoint, upperEndpoint;
  for(int i = 0; i < nrows; i ++) {
    for(int j = 0; j < ncols; j ++) {
      
      // standard error calculation
      standardError(i, j) = sigma(i, j) / sqrt(double(n));
      
      // marginal error calculation
      marginOfError = zscore * standardError(i, j);
      
      // lower and upper end point calculation
      lowerEndpoint = transMatr(i, j) - marginOfError;
      upperEndpoint = transMatr(i, j) + marginOfError;
      
      // taking care that upper and lower end point should be between 0(included) and 1(included)
      lowerEndpointMatr(i, j) = (lowerEndpoint > 1.0) ? 1.0 : ((0.0 > lowerEndpoint) ? 0.0 : lowerEndpoint);
      upperEndpointMatr(i, j) = (upperEndpoint > 1.0) ? 1.0 : ((0.0 > upperEndpoint) ? 0.0 : upperEndpoint);
    }
  }
  
  // set the rows and columns name as states names
  standardError.attr("dimnames") = upperEndpointMatr.attr("dimnames") 
              = lowerEndpointMatr.attr("dimnames") = transMatr.attr("dimnames"); 
  
  // return a list of important results
  List out = List::create(_["estimate"] = estimate,
                          _["standardError"] = standardError,
                          _["confidenceInterval"] = List::create(_["confidenceLevel"] = confidencelevel, 
                                                                 _["lowerEndpointMatrix"] = lowerEndpointMatr,
                                                                 _["upperEndpointMatrix"] = upperEndpointMatr),
                          _["bootStrapSamples"] = pmsBootStrapped
		                     ); 

  return out;
}

// convert matrix data to transition probability matrix
S4 _matr2Mc(CharacterMatrix matrData, double laplacian = 0) {
  
  // dimension of input matrix
  int nRows = matrData.nrow(), nCols = matrData.ncol();

  // set of states
  std::set<std::string> uniqueVals;
  
  // populate uniqueVals set
  for(int i = 0; i < nRows; i++) 
  	for(int j = 0; j < nCols; j++) 
		  uniqueVals.insert((std::string)matrData(i, j));	

  // unique states
  int usize = uniqueVals.size();
  
  // matrix of dimension usize
  NumericMatrix contingencyMatrix (usize);
  
  // state names as rows name and columns name
  contingencyMatrix.attr("dimnames") = List::create(uniqueVals, uniqueVals); 
  
  // iterator for set of states
  std::set<std::string>::iterator it;
  
  // populate contingency matrix
  int stateBegin = 0, stateEnd = 0;
  for(int i = 0; i < nRows; i ++) {
    for(int j = 1; j < nCols; j ++) {
      
      // row and column number of begin state and end state
      int k = 0;
      for(it = uniqueVals.begin(); it != uniqueVals.end(); ++it, k++) {
        if(*it == (std::string)matrData(i, j-1)) {
          stateBegin = k;
        }
        
        if(*it == (std::string)matrData(i,j)) {
          stateEnd = k;
        }
      }
      
      contingencyMatrix(stateBegin,stateEnd)++;
    }
  }

  // add laplacian correction if needed
  for(int i = 0; i < usize; i++) {
	  double rowSum = 0;
	  for(int j = 0; j < usize; j++) {
    		contingencyMatrix(i,j) += laplacian;
    		rowSum += contingencyMatrix(i, j);
  	}
  	
  	// get the transition matrix and a DTMC
	  for(int j = 0; j < usize; j ++) {
	    contingencyMatrix(i,j) /= rowSum;
	  }
    		
  }
  
  // markovchain object
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
        if(transMatr(i, j) < 0. || transMatr(i, j) > 1.)
          stop("The entries in the transition matrix must each belong to the interval [0, 1]");
        else
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

// main function for fitting markov chains
// [[Rcpp::export]]
List markovchainFit(SEXP data, String method = "mle", bool byrow = true, int nboot = 10, double laplacian = 0
            , String name = "", bool parallel = false, double confidencelevel = 0.95
            , NumericMatrix hyperparam = NumericMatrix()) {
  
  // list to store the output
  List out;
  
  // if input data is data frame or matrix
  if(Rf_inherits(data, "data.frame") || Rf_inherits(data, "matrix")) { 
    
    // store input data in mat
    CharacterMatrix mat;
    
    // if data is a data.frame force it to matrix
  	if(Rf_inherits(data, "data.frame")) {
  	  DataFrame df(data);
  	  
  	  // matrix : no of rows = no of rows in df : same for number of columns
  	  mat = CharacterMatrix(df.nrows(), df.size());
  	  for(int i = 0; i < df.size(); i++) {
  	    mat(_,i) = CharacterVector(df[i]);
  	  }
 	  } 
  	else {
		  mat = data;
	  }
  	
  	// byrow assumes distinct observations (trajectiories) are per row
  	// otherwise transpose
  	if(!byrow) {
  	  mat = _transpose(mat); 
  	}
  	
   	S4 outMc =_matr2Mc(mat, laplacian);
   	out = List::create(_["estimate"] = outMc);
  } 
  else {
    
    if(method == "mle") {
      out = _mcFitMle(data, byrow, confidencelevel); 
    }
    
    if(method == "bootstrap") {
      out = _mcFitBootStrap(data, nboot, byrow, parallel, confidencelevel);
    }
  
    if(method == "laplace") {
      out = _mcFitLaplacianSmooth(data, byrow, laplacian);
    }
    
    if(method == "map") {
      out = _mcFitMap(data, byrow, confidencelevel, hyperparam);
    }
  }
  
  // markovchain object
  S4 estimate = out["estimate"];
  if(name != "") {
    estimate.slot("name") = name; 
  }
  
  // transition matrix
  NumericMatrix transMatr = estimate.slot("transitionMatrix");
  
  // data is neither data frame nor matrix
  if(!Rf_inherits(data, "data.frame") && !Rf_inherits(data, "matrix")) 
    out["logLikelihood"] = _loglikelihood(data, transMatr);
  
  estimate.slot("states") = rownames(transMatr);
  out["estimate"] = estimate;
  
  return out;
}