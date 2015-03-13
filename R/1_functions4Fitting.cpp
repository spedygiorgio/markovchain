#include <Rcpp.h>
#include <set>
#include <algorithm>
#include <random>
#include <cstdlib>
#include <ctime>
#include <math.h>
using namespace Rcpp;

template <typename T>
T transpose(const T & m){      
    int k = m.rows(), n = m.cols();
    T z(n, k);
    int sz = n*k-1;
    typename T::const_iterator mit;
    typename T::iterator zit;
    for(mit = m.begin(), zit = z.begin(); mit != m.end(); mit++, zit += n){
        if (zit >= z.end()) zit -= sz;
        *zit = *mit;
    }
    return(z);
}

template <typename T> 
T mean(std::vector<T> v){
	T mu = 0.;
	for(int i = 0; i < v.size(); i++)
		mu += v[i];
	return mu / v.size();
}

template <typename T>
T sd(std::vector<T> v){
	T sd = 0., mu = mean(v);
	for(int i = 0; i < v.size(); i++)
		sd += (v[i] - mu) * (v[i] - mu);
	return sqrt(sd / (v.size() - 1));
}

NumericVector rowSumsC(NumericMatrix x){
	int nrow = x.nrow(), ncol = x.ncol();
	NumericVector out(nrow);

	for(int i = 0; i < nrow; i++){
		double total = 0;
		for(int j = 0; j < ncol; j++){
		  	total += x(i, j);
		}
		out[i] = total;
	}
	return out;
}

/*
* core function to get sequence matrix
*/

// [[Rcpp::export]]
NumericMatrix createSequenceMatrix(CharacterVector stringchar, bool toRowProbs = 0, bool sanitize = 1, double laplacian = 0.0){
	// extract markov chain state names and store them in a set
	std::set<std::string> uniqueVals;
	for(int i = 0; i < stringchar.size(); i++)	
		uniqueVals.insert(std::string(stringchar[i]));
	// now store all in a sorted fasion in a vector
	std::vector<std::string> keys;
	for(std::set<std::string>::iterator it = uniqueVals.begin(); it != uniqueVals.end(); it++)	
		keys.push_back(*it);

	int sizeMatr = keys.size();
	NumericMatrix freqMatrix = NumericMatrix(sizeMatr, sizeMatr);
	std::fill(freqMatrix.begin(), freqMatrix.end(), 0);
	freqMatrix.attr("rownames") = freqMatrix.attr("colnames") = List::create(keys);

	// populate the frequency matrix
	for(int i = 0; i < stringchar.size() - 1; i++){
		std::string startState = std::string(stringchar(i));
		std::string endState = std::string(stringchar(i+1));
		std::vector<std::string>::iterator whichRow = std::lower_bound(keys.begin(), keys.end(), startState);
		std::vector<std::string>::iterator whichCol = std::lower_bound(keys.begin(), keys.end(), endState);
		freqMatrix(whichRow - keys.begin(), whichCol - keys.begin())++;
	}

	NumericVector rowSum = rowSumsC(freqMatrix);

	// sanitize rows whose row sums are zeros
	if(sanitize)
		for(int i = 0; i < sizeMatr; i++)
			if(!rowSum[i])
				for(int j = 0; j < sizeMatr; j++)	
					freqMatrix(i, j) = 1 / sizeMatr;

	// laplacian smoothening here
	for(int i = 0; i < sizeMatr; i++)
		for(int j = 0; j < sizeMatr; j++)
			freqMatrix(i, j) += laplacian;

	rowSum = rowSumsC(freqMatrix);

	if(toRowProbs)
		for(int i = 0; i < sizeMatr; i++)
			for(int j = 0; j < sizeMatr; j++)
				freqMatrix(i, j) /= rowSum(i);
	
	return freqMatrix;
}

/*
* maximum likelihood mc fitter
*/

// [[Rcpp::export]]
List mcFitMle(CharacterVector stringchar, bool byrow){
	NumericMatrix initialMatrix = createSequenceMatrix(stringchar, 1, 1, 0.0);

	if(!byrow) initialMatrix = transpose(initialMatrix);
	S4 outMc("markovchain");
	outMc.slot("transitionMatrix") = initialMatrix;
	outMc.slot("name") = "MLE Fit";

	List out = List::create(Named("estimate") = outMc);
	return out;
}

/*
* maximum likelihood mc fitter with laplacian smoothening
*/

// [[Rcpp::export]]
List mcFitLaplacianSmooth(CharacterVector stringchar, bool byrow, double laplacian = 0.01){
	NumericMatrix initialMatrix = createSequenceMatrix(stringchar, 1, 1, laplacian);
	
	if(!byrow) initialMatrix = transpose(initialMatrix);
	S4 outMc("markovchain");
	outMc.slot("transitionMatrix") = initialMatrix;
	outMc.slot("name") = "MLE Fit";

	List out = List::create(Named("estimate") = outMc);
	return out;
}

/*
* function which generates the bootstrap character sequences
*/

// [[Rcpp::export]]
List bootstrapCharacterSequences(CharacterVector stringchar, int n){
	// n is the number of bootstrap sets (nboot)
	// seed the random number generator
	srand(time(0));
	std::default_random_engine generator;
  	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	NumericMatrix contingencyMatrix = createSequenceMatrix(stringchar, 1);
	List samples(n);
	List _rownames = contingencyMatrix.attr("rownames");
	CharacterVector rownames = as<CharacterVector>(_rownames(0));
	int size = stringchar.size();

	for(int i = 0; i < n; i++){
		std::vector<std::string> charseq;
		int idx = rand() % contingencyMatrix.ncol();
		std::string item = std::string(rownames(idx));
		charseq.push_back(item);

		for(int j = 1; j < size; j++){
			NumericVector probVector = contingencyMatrix(idx, _);
			std::vector<double> probs;
			probs.push_back(probVector(0));
			for(int k = 1; k < contingencyMatrix.ncol(); k++)	
				probs.push_back(probs[k-1] + probVector(k));

			double res = distribution(generator);
			charseq.push_back(std::string(rownames(int(std::lower_bound(probs.begin(), probs.end(), res) - probs.begin()))));
		}
		samples(i) = charseq;
	}
	return samples;
}

/*
* from the bootstrap data set, estimate the parameters
*/

// [[Rcpp::export]]
List fromBoot2Estimate(List listMatr){
	int sampleSize = listMatr.size();
	std::vector<NumericMatrix> arrayMatr;
	for(int i = 0; i < sampleSize; i++)
		arrayMatr.push_back(as<NumericMatrix>(listMatr(i)));
	int matrDim = arrayMatr[0].nrow();
	
	NumericMatrix matrMean = NumericMatrix(matrDim, matrDim);
	NumericMatrix matrSd = NumericMatrix(matrDim, matrDim);

	for(int i = 0; i < matrDim; i++){
		for(int j = 0; j < matrDim; j++){
			std::vector<double> probsEst;
			for(int k = 0; k < sampleSize; k++)
				probsEst.push_back(arrayMatr[k](i, j));
			matrMean(i, j) = mean(probsEst);
			matrSd(i, j) = sd(probsEst);
		}
	}

	return List::create(Named("estMu") = matrMean, Named("estSigma") = matrSd);
}

// [[Rcpp::export]]
S4 matr2Mc(CharacterMatrix matrData, double laplacian = 0){
	int nCols = matrData.ncol(), nRows = matrData.nrow();

	// store all the unique values of the matrix in a set (stored in sorted order)
	std::set<std::string> uniqueVals;
	for(int i = 0; i < nRows; i++)
		for(int j = 0; j < nCols; j++)
			uniqueVals.insert(std::string(matrData(i, j)));

	NumericMatrix contingencyMatrix(uniqueVals.size(), uniqueVals.size());
	std::fill(contingencyMatrix.begin(), contingencyMatrix.end(), 0);

	// extract all sorted unique values (keys)
	std::vector<std::string> keys;
	for(std::set<std::string>::iterator it = uniqueVals.begin(); it != uniqueVals.end(); it++)	keys.push_back(*it);
	contingencyMatrix.attr("rownames") = contingencyMatrix.attr("colnames") = List::create(keys);

	// populate contingency matrix
	for(int i = 0; i < nRows; i++){
		for(int j = 1; j < nCols; j++){
			std::string startState = std::string(matrData(i, j-1));
			std::string endState = std::string(matrData(i, j));
			std::vector<std::string>::iterator whichRow = std::lower_bound(keys.begin(), keys.end(), startState);
			std::vector<std::string>::iterator whichCol = std::lower_bound(keys.begin(), keys.end(), endState);
			contingencyMatrix(whichRow - keys.begin(), whichCol - keys.begin())++;
		}
	}

	// add laplacian
	for(int i = 0; i < nRows; i++)
		for(int j = 0; j < nCols; j++)
			contingencyMatrix(i, j) += laplacian;
	// compute transition matrix
	NumericVector rowSum = rowSumsC(contingencyMatrix);
	NumericMatrix transitionMatrix(uniqueVals.size(), uniqueVals.size());
	std::fill(transitionMatrix.begin(), transitionMatrix.end(), 0);
	for(int i = 0; i < nRows; i++)
		for(int j = 0; j < nCols; j++)
			transitionMatrix(i, j) = contingencyMatrix(i, j) / rowSum(i);

	S4 outMc("markovchain");
	outMc.slot("transitionMatrix") = transitionMatrix;

	return outMc;
}