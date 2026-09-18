#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <armadillo>
#include <Rcpp.h>
#include <cmath>
#include <limits>
using namespace Rcpp;
using namespace RcppArmadillo;
using namespace arma;
using namespace std;



// [[Rcpp::export(.ExpectedTimeRCpp)]]
NumericVector ExpectedTimeRcpp(NumericMatrix x, NumericVector y) {
  NumericVector out;
  int size = x.nrow();
  arma::mat T = arma::zeros(size, size);

  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      T(i, j) = x(i, j);
    
  arma::vec c = arma::zeros(size);
  
  for (int i = 0; i < size; i++)
    c[i] = y[i];

  out = wrap(solve(T,c));

  return out;
}


// [[Rcpp::export(.probabilityatTRCpp)]]
NumericMatrix probabilityatTRCpp(NumericMatrix y) {
  
  int size = y.nrow();
  NumericMatrix out(size,size);

  arma::mat T = arma::zeros(size, size);
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      T(i, j) = y(i, j);

  T = expmat(T);

  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      out(i, j) = T(i, j);

  return out;
}



// [[Rcpp::export(.impreciseProbabilityatTRCpp)]]
NumericVector impreciseProbabilityatTRCpp(S4 C, int i,int t, int s, double error) {
  if (!std::isfinite(error) || error <= 0.0)
    stop("error must be finite and positive");
  if (s <= t) stop("s must be greater than t");
  CharacterVector states = C.slot("states");
  int noOfstates = states.size();
  
  NumericMatrix non_Q = C.slot("Q");
  NumericMatrix non_range = C.slot("range");
  
  arma::mat Q = arma::zeros(noOfstates, noOfstates);
  arma::mat range = arma::zeros(noOfstates, 2);
  
  // initialises armadillo matrices
  for (int p = 0; p < noOfstates; p++)
    for (int q = 0; q < noOfstates; q++)
      Q(p, q) = non_Q(p, q);
  
  // initialises armadillo matrices
  for (int p = 0; p < noOfstates; p++)
    for (int q = 0; q < 2; q++)
      range(p, q) = non_range(p, q);
  
  // initialses value of norm of Q
  double QNorm = -1.0;
  
  // calculates norm of Q
  for (int p =0; p < noOfstates; p++) {
    float sum = 0.0;
    
    for (int q = 0; q < noOfstates; q++) {
      if (Q(p, q) >= 0)
        sum = sum + Q(p, q);
      else
        sum = sum + -Q(p, q);
    }
    
    if (sum * range(p, 1) > QNorm)
      QNorm = sum * range(p, 1);
  }
  
  // Calculate the number of iterations in floating point, validate it before
  // conversion to int, and guarantee at least one iteration.
  const double horizon = static_cast<double>(s - t);
  const double candidate = std::max(
    horizon * QNorm,
    horizon * horizon * QNorm * QNorm / (2.0 * error));
  if (!std::isfinite(candidate) ||
      candidate > static_cast<double>(std::numeric_limits<int>::max()))
    stop("requested accuracy requires too many iterations");
  const int n = std::max(1, static_cast<int>(std::ceil(candidate)));
  
  // sets delta value
  const double delta = horizon / static_cast<double>(n);
  
  // declares and initialises initial f
  arma::vec Ii(noOfstates);
  
  for (int p = 0; p < noOfstates; p++)
    Ii[p] = 0;
  
  Ii[i-1] = 1;
  
  
  // calculation of Qgx vector
  arma::vec values =  Q * Ii;
  arma::vec Qgx(noOfstates);
  
  for (int p = 0; p < noOfstates; p++)
    Qgx[p] = 0;
  
  for (int p = 0; p < noOfstates; p++) {
    if (values[p] * range(p, 0) < values[p] * range(p, 1))
      Qgx[p] = values[p] * range(p,0);
    else
      Qgx[p] = values[p] * range(p,1);
  }
  
  Qgx = delta * Qgx;
  
  // runs n-1 iterations according to the algorithm
  // Qgx_i = Qgx_{i-1} + delta*Q*Qgx_{i-1}
  Qgx = Qgx + Ii;
  
  for (int iter = 0; iter < n - 1; iter++) {
    arma::vec temp = Qgx;
    values = Q * Qgx;
    
    for (int p = 0; p < noOfstates; p++) {
      // calculating keeping in mind the lower opertaotr values
      if (values[p] * range(p,0) < values[p] * range(p,1))
        Qgx[p] = values[p] * range(p,0);
      else
        Qgx[p] = values[p] * range(p,1);
    }
    
    Qgx = delta * Qgx;
    Qgx = temp + Qgx;
  }
  
  NumericVector out;
  
  for (int p = 0; p < noOfstates; p++)
    out.push_back(Qgx[p]);
  
  return out;
}


