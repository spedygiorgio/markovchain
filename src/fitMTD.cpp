#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// EM algorithm for the mixture transition distribution (MTD) model of
// Raftery (1985) with non-negative lag weights (Lebre and Bourguignon, 2008):
//
//   P(X_t = j | X_{t-1} = i_1, ..., X_{t-k} = i_k) = sum_g lambda_g Q(i_g, j),
//
// with a single row-stochastic matrix Q (rows are the departure state) shared
// by all lags. The latent variable is the lag that generated X_t; the E-step
// computes its posterior probability and the M-step updates lambda as the
// average posterior of each lag and Q as posterior-weighted transition counts.
// Each iteration does not decrease the log-likelihood.
//
// patterns: m x (k + 1) matrix of 0-based state indices, column 0 is X_t and
//           column g is X_{t-g}; each row is a distinct observed pattern
// counts:   number of times each pattern occurs
// lambda0, Q0: starting values (Q0 rows sum to one)
// The log-likelihood returned is the one of the parameters returned.
// [[Rcpp::export(.mtdEM)]]
List mtdEM(IntegerMatrix patterns, NumericVector counts, NumericVector lambda0,
           NumericMatrix Q0, double tol, int maxit) {
  const int m = patterns.nrow();
  const int k = patterns.ncol() - 1;
  const int r = Q0.nrow();
  if (k < 1 || lambda0.size() != k || Q0.ncol() != r || counts.size() != m)
    stop("inconsistent dimensions");

  NumericVector lambda = clone(lambda0);
  NumericMatrix Q = clone(Q0);
  NumericVector lambdaNew(k);
  NumericMatrix Qnew(r, r);
  std::vector<double> w(k);
  double total = 0.0;
  for (int p = 0; p < m; ++p) total += counts[p];

  double logLik = R_NegInf, logLikOld = R_NegInf;
  bool converged = false;
  int iter = 0;
  for (iter = 1; iter <= maxit; ++iter) {
    // E-step and sufficient statistics at the current parameters
    std::fill(lambdaNew.begin(), lambdaNew.end(), 0.0);
    std::fill(Qnew.begin(), Qnew.end(), 0.0);
    logLik = 0.0;
    for (int p = 0; p < m; ++p) {
      const int to = patterns(p, 0);
      double mix = 0.0;
      for (int g = 0; g < k; ++g) {
        w[g] = lambda[g] * Q(patterns(p, g + 1), to);
        mix += w[g];
      }
      if (!(mix > 0.0)) {
        logLik = R_NegInf;
        break;
      }
      logLik += counts[p] * std::log(mix);
      for (int g = 0; g < k; ++g) {
        const double post = counts[p] * w[g] / mix;
        lambdaNew[g] += post;
        Qnew(patterns(p, g + 1), to) += post;
      }
    }
    if (!std::isfinite(logLik))
      stop("an observed pattern has probability zero at the starting values");
    if (logLik - logLikOld <= tol * (std::fabs(logLik) + tol)) {
      converged = true;
      break;
    }
    if (iter == maxit) break;
    logLikOld = logLik;
    // M-step
    for (int g = 0; g < k; ++g) lambda[g] = lambdaNew[g] / total;
    for (int i = 0; i < r; ++i) {
      double rowTotal = 0.0;
      for (int j = 0; j < r; ++j) rowTotal += Qnew(i, j);
      // a state that never appears in a conditioning position with positive
      // weight keeps its row: it does not enter the likelihood
      if (rowTotal > 0.0)
        for (int j = 0; j < r; ++j) Q(i, j) = Qnew(i, j) / rowTotal;
    }
  }
  return List::create(_["lambda"] = lambda, _["Q"] = Q, _["logLik"] = logLik,
                      _["iterations"] = iter > maxit ? maxit : iter,
                      _["converged"] = converged);
}
