#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// This file implements the EM algorithm (Expectation-Maximisation) of the MTD model in C++ for speed;
// R/fitMTD.R calls it. EM in two lines: we do not know which lag 'generated' each observation
// (latent variable). E-step: compute the posterior probability of each lag for every observation.
// M-step: update the parameters as weighted averages using those probabilities. Repeat until the
// log-likelihood stops increasing (it never decreases).
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
// mtdEM(): iterate E-step and M-step until convergence (relative change of the log-likelihood below tol)
// or maxit iterations. The data come as distinct patterns with counts to avoid looping over
// identical observations.
// [[Rcpp::export(.mtdEM)]]
List mtdEM(IntegerMatrix patterns, NumericVector counts, NumericVector lambda0,
           NumericMatrix Q0, double tol, int maxit) {
  const int m = patterns.nrow();
  const int k = patterns.ncol() - 1;
  const int r = Q0.nrow();
  if (k < 1 || lambda0.size() != k || Q0.ncol() != r || counts.size() != m)
    stop("inconsistent dimensions");

  // Work on copies (clone): Rcpp vectors share memory with the R object, and we must not modify the input.
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
    // E-step. For each pattern: w[g] = lambda[g] * Q(previous state at lag g, current state) is the joint
    // probability that lag g generated the observation; mix = sum over g = probability of the pattern.
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
      // Log-likelihood = sum counts * log(mix). Posterior of lag g = w[g] / mix.
      logLik += counts[p] * std::log(mix);
      for (int g = 0; g < k; ++g) {
        const double post = counts[p] * w[g] / mix;
        lambdaNew[g] += post;
        Qnew(patterns(p, g + 1), to) += post;
      }
    }
    if (!std::isfinite(logLik))
      stop("an observed pattern has probability zero at the starting values");
    // Convergence test on the relative increase of the log-likelihood.
    if (logLik - logLikOld <= tol * (std::fabs(logLik) + tol)) {
      converged = true;
      break;
    }
    if (iter == maxit) break;
    logLikOld = logLik;
    // M-step: lambda = average posterior; Q(i, .) = posterior-weighted transition counts, normalised by row.
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

// mtdStep(): ONE iteration of the same algorithm, exposed to R for the SQUAREM acceleration
// (see .accelerateEM in R/fitMTD.R). It also returns the log-likelihood of the INPUT parameters.
// ---------------------------------------------------------------------------
// One EM step for the MTD model with a single matrix shared by all lags, used
// by the SQUAREM acceleration written in R. It returns F(theta), the
// parameters after one EM iteration, and the log-likelihood at the INPUT
// parameters (-Inf, without error, if an observed pattern has probability
// zero, so that a rejected extrapolation can be backtracked). The updates are
// those of mtdEM().
// [[Rcpp::export(.mtdStep)]]
List mtdStep(IntegerMatrix patterns, NumericVector counts, NumericVector lambda,
             NumericMatrix Q) {
  const int m = patterns.nrow();
  const int k = patterns.ncol() - 1;
  const int r = Q.nrow();
  if (k < 1 || lambda.size() != k || Q.ncol() != r || counts.size() != m)
    stop("inconsistent dimensions");
  NumericVector lambdaNew(k);
  NumericMatrix Qnew(r, r);
  std::vector<double> w(k);
  double total = 0.0, logLik = 0.0;
  for (int p = 0; p < m; ++p) total += counts[p];
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
  NumericVector lambdaOut(k);
  NumericMatrix Qout = clone(Q);
  if (std::isfinite(logLik)) {
    for (int g = 0; g < k; ++g) lambdaOut[g] = lambdaNew[g] / total;
    for (int i = 0; i < r; ++i) {
      double rowTotal = 0.0;
      for (int j = 0; j < r; ++j) rowTotal += Qnew(i, j);
      if (rowTotal > 0.0)
        for (int j = 0; j < r; ++j) Qout(i, j) = Qnew(i, j) / rowTotal;
    }
  } else {
    lambdaOut = clone(lambda);
  }
  return List::create(_["lambda"] = lambdaOut, _["Q"] = Qout, _["logLik"] = logLik);
}

// MTDg model: same EM but each lag g has its own matrix Q_g. The three helpers below perform the E-step,
// the M-step and the copying of the matrices; mtdgEM and mtdgStep use them.
// ---------------------------------------------------------------------------
// EM for the MTDg model, in which every lag g has its own row-stochastic
// matrix Q_g (rows are the departure state):
//
//   P(X_t = j | X_{t-1} = i_1, ..., X_{t-k} = i_k) = sum_g lambda_g Q_g(i_g, j).
//
// The latent variable is again the lag that generated X_t. The E-step gives its
// posterior probability; the M-step updates lambda_g as the average posterior of
// lag g and each Q_g as the posterior-weighted transition counts of its own lag.
// Every iteration does not decrease the log-likelihood.
// Arguments as in mtdEM(), with the list Qlist0 of the k starting matrices.

// E-step of MTDg: accumulates the weights of the lags and the posterior-weighted counts for each lag's matrix.
static double mtdgEStep(const IntegerMatrix& patterns, const NumericVector& counts,
                        const NumericVector& lambda, const std::vector<NumericMatrix>& Q,
                        NumericVector& lambdaNew, std::vector<NumericMatrix>& Qnew,
                        std::vector<double>& w) {
  const int m = patterns.nrow();
  const int k = lambda.size();
  std::fill(lambdaNew.begin(), lambdaNew.end(), 0.0);
  for (int g = 0; g < k; ++g) std::fill(Qnew[g].begin(), Qnew[g].end(), 0.0);
  double logLik = 0.0;
  for (int p = 0; p < m; ++p) {
    const int to = patterns(p, 0);
    double mix = 0.0;
    for (int g = 0; g < k; ++g) {
      w[g] = lambda[g] * Q[g](patterns(p, g + 1), to);
      mix += w[g];
    }
    if (!(mix > 0.0)) return R_NegInf;
    logLik += counts[p] * std::log(mix);
    for (int g = 0; g < k; ++g) {
      const double post = counts[p] * w[g] / mix;
      lambdaNew[g] += post;
      Qnew[g](patterns(p, g + 1), to) += post;
    }
  }
  return logLik;
}

// M-step of MTDg: normalise.
static void mtdgMStep(const NumericVector& lambdaNew, const std::vector<NumericMatrix>& Qnew,
                      double total, NumericVector& lambda, std::vector<NumericMatrix>& Q) {
  const int k = lambda.size();
  const int r = Q[0].nrow();
  for (int g = 0; g < k; ++g) {
    lambda[g] = lambdaNew[g] / total;
    for (int i = 0; i < r; ++i) {
      double rowTotal = 0.0;
      for (int j = 0; j < r; ++j) rowTotal += Qnew[g](i, j);
      // a state never observed at this lag keeps its row: it does not enter the likelihood
      if (rowTotal > 0.0)
        for (int j = 0; j < r; ++j) Q[g](i, j) = Qnew[g](i, j) / rowTotal;
    }
  }
}

// Copy a list of R matrices into C++ vectors of matrices, checking the dimensions.
static std::vector<NumericMatrix> mtdgCopyMatrices(const List& Qlist, int k, int r) {
  if (Qlist.size() != k) stop("inconsistent dimensions");
  std::vector<NumericMatrix> Q(k);
  for (int g = 0; g < k; ++g) {
    NumericMatrix M = as<NumericMatrix>(Qlist[g]);
    if (M.nrow() != r || M.ncol() != r) stop("inconsistent dimensions");
    Q[g] = clone(M);
  }
  return Q;
}

static List mtdgWrapMatrices(const std::vector<NumericMatrix>& Q) {
  List out(Q.size());
  for (size_t g = 0; g < Q.size(); ++g) out[g] = Q[g];
  return out;
}

// mtdgEM(): iterate until convergence, as mtdEM().
// [[Rcpp::export(.mtdgEM)]]
List mtdgEM(IntegerMatrix patterns, NumericVector counts, NumericVector lambda0,
            List Qlist0, double tol, int maxit) {
  const int m = patterns.nrow();
  const int k = patterns.ncol() - 1;
  if (k < 1 || lambda0.size() != k || counts.size() != m || Qlist0.size() != k)
    stop("inconsistent dimensions");
  const int r = as<NumericMatrix>(Qlist0[0]).nrow();
  std::vector<NumericMatrix> Q = mtdgCopyMatrices(Qlist0, k, r);
  std::vector<NumericMatrix> Qnew(k);
  for (int g = 0; g < k; ++g) Qnew[g] = NumericMatrix(r, r);
  NumericVector lambda = clone(lambda0);
  NumericVector lambdaNew(k);
  std::vector<double> w(k);
  double total = 0.0;
  for (int p = 0; p < m; ++p) total += counts[p];

  double logLik = R_NegInf, logLikOld = R_NegInf;
  bool converged = false;
  int iter = 0;
  for (iter = 1; iter <= maxit; ++iter) {
    logLik = mtdgEStep(patterns, counts, lambda, Q, lambdaNew, Qnew, w);
    if (!std::isfinite(logLik))
      stop("an observed pattern has probability zero at the starting values");
    if (logLik - logLikOld <= tol * (std::fabs(logLik) + tol)) {
      converged = true;
      break;
    }
    if (iter == maxit) break;
    logLikOld = logLik;
    mtdgMStep(lambdaNew, Qnew, total, lambda, Q);
  }
  return List::create(_["lambda"] = lambda, _["Q"] = mtdgWrapMatrices(Q),
                      _["logLik"] = logLik,
                      _["iterations"] = iter > maxit ? maxit : iter,
                      _["converged"] = converged);
}

// One EM step of the MTDg model: F(theta) and the log-likelihood at the input
// parameters (-Inf if an observed pattern has probability zero).
// [[Rcpp::export(.mtdgStep)]]
List mtdgStep(IntegerMatrix patterns, NumericVector counts, NumericVector lambda,
              List Qlist) {
  const int m = patterns.nrow();
  const int k = patterns.ncol() - 1;
  if (k < 1 || lambda.size() != k || counts.size() != m || Qlist.size() != k)
    stop("inconsistent dimensions");
  const int r = as<NumericMatrix>(Qlist[0]).nrow();
  std::vector<NumericMatrix> Q = mtdgCopyMatrices(Qlist, k, r);
  std::vector<NumericMatrix> Qnew(k);
  for (int g = 0; g < k; ++g) Qnew[g] = NumericMatrix(r, r);
  NumericVector lambdaNew(k);
  NumericVector lambdaOut = clone(lambda);
  std::vector<double> w(k);
  double total = 0.0;
  for (int p = 0; p < m; ++p) total += counts[p];
  const double logLik = mtdgEStep(patterns, counts, lambda, Q, lambdaNew, Qnew, w);
  if (std::isfinite(logLik)) mtdgMStep(lambdaNew, Qnew, total, lambdaOut, Q);
  return List::create(_["lambda"] = lambdaOut, _["Q"] = mtdgWrapMatrices(Q),
                      _["logLik"] = logLik);
}
// The code that follows (clean_nas and _mcFitMap) belongs to mapFitFunctions.h, see the comments there.
