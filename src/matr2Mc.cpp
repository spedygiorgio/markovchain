#include <Rcpp.h>
#include <iterator>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector rowSumsRcpp(NumericMatrix matr) {
	int nrow = matr.nrow(), ncol = matr.ncol();
	NumericVector rowSums(nrow);
	for (int i = 0; i < nrow; ++i){
		double sum = 0;
		for (int j = 0; j < ncol; ++j){
			sum += matr(i, j);
		}
		rowSums[i] = sum;
	}

	return rowSums;
}
// [[Rcpp::export]]
NumericMatrix matr2Mc(CharacterMatrix matrData, double laplacian = 0) {
	//find unique values scanning the matrix
	std::set<std::string> uniqueVals;
	int nrow = matrData.nrow(), ncol = matrData.ncol();

	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			uniqueVals.insert((std::string)matrData(i, j));
		}
	}
	//create a contingency matrix
	NumericMatrix contingencyMatrix(uniqueVals.size(), uniqueVals.size());
	std::fill(contingencyMatrix.begin(), contingencyMatrix.end(), 0);

	//fill the contingency matrix
	// Rcpp::List dimnms4c = Rcpp::List::create(uniqueVals, uniqueVals);
	// contingencyMatrix.attr("dimnames") = dimnms4c;
	for (int i = 0; i < nrow; ++i) {
		for (int j = 1; j < ncol; ++j) {
			std::set<std::string>::iterator itForRow;
			itForRow = std::find (uniqueVals.begin(), uniqueVals.end(), (std::string)matrData(i, j - 1));
			int whichRow = std::distance(uniqueVals.begin(), itForRow);
			std::set<std::string>::iterator itForCols;
			itForCols = find (uniqueVals.begin(), uniqueVals.end(), (std::string)matrData(i, j));
			int whichCols = std::distance(uniqueVals.begin(), itForCols);
			contingencyMatrix(whichRow, whichCols) = contingencyMatrix(whichRow, whichCols) + 1;
		}
	}


	//add laplacian correction if needed
	if (laplacian != 0) {
		for (int i = 0; i < nrow; ++i) {
			for (int j = 0; j < ncol; ++j) {
				contingencyMatrix(i, j) += laplacian;
			}
		}
	}

	//get a transition matrix and a DTMC
	NumericMatrix transitionMatrix(uniqueVals.size(), uniqueVals.size());
	std::fill(transitionMatrix.begin(), transitionMatrix.end(), 0);
	Rcpp::List dimnms4t = Rcpp::List::create(uniqueVals, uniqueVals);
	transitionMatrix.attr("dimnames") = dimnms4t;
	NumericVector rowSums = rowSumsRcpp(contingencyMatrix);
	for (int i = 0; i < nrow; ++i) {
		for (int j = 0; j < ncol; ++j) {
			transitionMatrix(i, j) = contingencyMatrix(i, j) / rowSums(i);
		}
	}



	return transitionMatrix;

}