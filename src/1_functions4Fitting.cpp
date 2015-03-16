#include <Rcpp.h>
#include <omp.h>
#include <unistd.h>

using namespace Rcpp;

template <typename T>
T _transpose(T & m) {      
  int k = m.rows(), n = m.cols();
  T z(n, k);
  z.attr("dimnames") = List::create(colnames(m), rownames(m)); 
  int sz1 = n*k-1;
  typename T::iterator mit, zit;
  for (mit = m.begin(), zit = z.begin(); mit != m.end(); mit++, zit += n) {
  	if (zit >= z.end()) zit -= sz1;
        *zit = *mit;
  }
  return (z);
}

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
NumericMatrix createSequenceMatrix_cpp(CharacterVector stringchar, bool toRowProbs=false, bool sanitize=true) {
  CharacterVector elements = unique(stringchar).sort();
  int sizeMatr = elements.size();
  
  NumericMatrix freqMatrix(sizeMatr);
  freqMatrix.attr("dimnames") = List::create(elements, elements); 
  CharacterVector rnames = rownames(freqMatrix);

  int posFrom, posTo;
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

List _mcFitMle(CharacterVector stringchar, bool byrow, double confidencelevel=95.0) {
  // get initialMatr and freqMatr 
  CharacterVector elements = unique(stringchar).sort();
  int sizeMatr = elements.size();
  
  NumericMatrix initialMatr(sizeMatr);
  NumericMatrix freqMatr(sizeMatr);
  initialMatr.attr("dimnames") = List::create(elements, elements); 

  int posFrom, posTo;
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
  double zscore = 1.96;
  if(confidencelevel == 99.9) zscore = 3.3;
  else if(confidencelevel == 99.0) zscore = 2.577;
  else if(confidencelevel == 98.5) zscore = 2.43;
  else if(confidencelevel == 97.5) zscore = 2.243;
  else if(confidencelevel == 90.0) zscore = 1.645;
  else if(confidencelevel == 85.0) zscore = 1.439;
  else if(confidencelevel == 75.0) zscore = 1.151;

  int n = stringchar.size();
  NumericMatrix lowerEndpointMatr = NumericMatrix(initialMatr.nrow(), initialMatr.ncol());
  NumericMatrix upperEndpointMatr = NumericMatrix(initialMatr.nrow(), initialMatr.ncol());

  double marginOfError, lowerEndpoint, upperEndpoint;
  for(int i = 0; i < initialMatr.nrow(); i ++) {
	for(int j = 0; j < initialMatr.ncol(); j ++) {
		marginOfError = zscore * initialMatr(i, j) / sqrt(freqMatr(i, j));
		lowerEndpoint = initialMatr(i, j) - marginOfError;
		upperEndpoint = initialMatr(i, j) + marginOfError;
		lowerEndpointMatr(i,j) = (lowerEndpoint > 1.0) ? 1.0 : ((0.0 > lowerEndpoint) ? 0.0 : lowerEndpoint);
		upperEndpointMatr(i,j) = (upperEndpoint > 1.0) ? 1.0 : ((0.0 > upperEndpoint) ? 0.0 : upperEndpoint);
  	}
  }
  lowerEndpointMatr.attr("dimnames") = List::create(elements, elements); 
  upperEndpointMatr.attr("dimnames") = List::create(elements, elements); 

  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = initialMatr;
  outMc.slot("name") = "MLE Fit";  
  
  return List::create(_["estimate"] = outMc
		, _["confidenceInterval"] = List::create(_["confidenceLevel"]=confidencelevel, 
							_["lowerEndpointMatrix"]=lowerEndpointMatr, 
							_["upperEndpointMatrix"]=upperEndpointMatr)
	);
}

List _mcFitLaplacianSmooth(CharacterVector stringchar, bool byrow, double laplacian=0.01) {
  NumericMatrix origNum = createSequenceMatrix_cpp(stringchar, false);
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
  NumericMatrix contingencyMatrix = createSequenceMatrix_cpp(stringchar);

  List samples, res;
  CharacterVector itemset = rownames(contingencyMatrix);
  int itemsetsize = itemset.size();

  Function sample("sample");
  srand(time(NULL));
  for(int i = 0; i < n; i ++) {
	CharacterVector charseq, resvec;	
	int rnd = rand()%itemsetsize;
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
		resvec = res[0]; ch = resvec[0];
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
  NumericMatrix matrMean(matrDim);
  NumericMatrix matrSd(matrDim);

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

List _mcFitBootStrap(CharacterVector data, int nboot=10, bool byrow=true, bool parallel=false) {
  List theList = _bootstrapCharacterSequences(data, nboot);
  int n = theList.size();
  List pmsBootStrapped(n);

  if(!parallel) { 
	for(int i = 0; i < n; i++) 
		pmsBootStrapped[i] = createSequenceMatrix_cpp(theList[i], true, true);
  } else {
	int cores = sysconf(_SC_NPROCESSORS_ONLN);
	omp_set_num_threads(cores);
	#pragma omp parallel for
	for(int i = 0; i < n; i ++) 
		pmsBootStrapped[i] = createSequenceMatrix_cpp(theList[i], true, true);
  }
  List estimateList = _fromBoot2Estimate(pmsBootStrapped);
  NumericMatrix transMatr = _toRowProbs(estimateList["estMu"]);

  S4 estimate("markovchain");
  estimate.slot("transitionMatrix") = transMatr;
  estimate.slot("byrow") = byrow;
  estimate.slot("name") = "BootStrap Estimate";  

  return List::create(_["estimate"] = estimate
		, _["standardError"] = estimateList["estSigma"]
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
  int stateBegin, stateEnd;
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
List markovchainFit_cpp(SEXP data, String method="mle", bool byrow=true, int nboot=10, double laplacian=0, String name="", bool parallel=false, double confidencelevel=0.95) {
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
    if(method == "mle") out = _mcFitMle(data, byrow);
    if(method == "bootstrap") out = _mcFitBootStrap(data, nboot, byrow, parallel);
    if(method == "laplace") out = _mcFitLaplacianSmooth(data, byrow, laplacian);
  }

  if(name != "") {
    S4 estimate = out["estimate"];
    estimate.slot("name") = name;
    out["estimate"] = estimate;
  }
  return out;
}

