#ifndef HELPERS_H
#define HELPERS_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;
using namespace std;

template <typename T>
T sortByDimNames(const T m);

double lbeta(double p, double q);

template <typename T>
T transposeMatrix(T & mat);

// BETAIN computes the incomplete Beta function ratio
double betain(double x, double p, double q, double beta);

// XINBTA computes inverse of the incomplete Beta function.
double xinbta(double p, double q, double beta, double alpha);

#endif
