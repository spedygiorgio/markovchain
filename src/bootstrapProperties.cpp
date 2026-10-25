// bootstrapProperties.cpp: the computational core of bootstrapProperties()
// and of the simulation in plotSequence() (R/bootstrapProperties.R and
// R/plotSequence.R).
//
// For a junior reader. A bootstrap replicate has two parts: drawing new
// transition counts and computing statistics of the chain they estimate. In
// R both are cheap per se, but every replicate also built a markovchain S4
// object (with its validity checks) and each of slem(), spectralGap() and
// relaxationTime() recomputed the same eigenvalues. Here the counts are
// drawn and the four built-in statistics are computed directly on plain
// matrices, with one eigen-decomposition per replicate.
//
// Random numbers come from R's generator (R::unif_rand, R::rmultinom), so
// set.seed() controls the results. The multinomial draws are made in the
// same order and with the same normalization as stats::rmultinom() in the
// R code they replace, so that bootstrap gives the same replicates as
// before for the same seed.
//
// The statistics reproduce slem(), spectralGap(), relaxationTime() and
// mixingTime() (R/spectralDiagnostics.R, R/dynamics.R,
// R/mixingDiagnostics.R): same definitions, same tolerances, and NA where
// the R function stops with an error (a reducible chain, or a periodic one
// for the mixing time). The package's tests compare the two on many chains.
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <queue>

using namespace Rcpp;

// Defined in probabilistic.cpp.
int periodKernel(NumericMatrix P);
arma::vec steadyStateErgodicMatrix(const arma::mat& submatrix);

namespace {

// Cumulative sums of the rows of a row-stochastic matrix, used to draw the
// next state by inversion: the next state is the first j with u < cum(i, j).
arma::mat rowCumulative(const arma::mat& P) {
  arma::mat cum = P;
  for (arma::uword i = 0; i < cum.n_rows; ++i) {
    for (arma::uword j = 1; j < cum.n_cols; ++j) {
      cum(i, j) += cum(i, j - 1);
    }
  }
  return cum;
}

// Draw a state from the row `row` of `cum` (or from a cumulative vector).
// A state of probability zero is never drawn: its cumulative value equals
// the previous one, so a smaller index is always found first. If rounding
// leaves the last cumulative value below u, the last state of positive
// probability is returned.
int drawFromCumulative(const double* cum, int n, int stride) {
  double u = R::unif_rand();
  for (int j = 0; j < n; ++j) {
    if (u < cum[j * stride]) {
      return j;
    }
  }
  for (int j = n - 1; j > 0; --j) {
    if (cum[j * stride] > cum[(j - 1) * stride]) {
      return j;
    }
  }
  return 0;
}

// Is the directed graph of the positive entries strongly connected? Every
// state must be reachable from state 0 and must reach it: one search on the
// graph and one on the reversed graph.
bool isStronglyConnected(const arma::mat& P) {
  const arma::uword n = P.n_rows;
  for (int reversed = 0; reversed < 2; ++reversed) {
    std::vector<bool> seen(n, false);
    std::queue<arma::uword> todo;
    seen[0] = true;
    todo.push(0);
    arma::uword count = 1;
    while (!todo.empty()) {
      arma::uword v = todo.front();
      todo.pop();
      for (arma::uword w = 0; w < n; ++w) {
        double edge = reversed ? P(w, v) : P(v, w);
        if (edge > 0 && !seen[w]) {
          seen[w] = true;
          ++count;
          todo.push(w);
        }
      }
    }
    if (count != n) {
      return false;
    }
  }
  return true;
}

// The four built-in statistics of a row-stochastic matrix, as computed by
// the R functions: which[k] is 0 for slem, 1 for spectralGap, 2 for
// relaxationTime and 3 for mixingTime; results go to out[k].
void spectralStatistics(const arma::mat& P, const IntegerVector& which,
                        double epsilon, int maxIter, double* out) {
  const int n = P.n_rows;
  const int k = which.size();
  // Everything is NA unless the statistic can be computed.
  for (int s = 0; s < k; ++s) {
    out[s] = NA_REAL;
  }
  // All four need an irreducible chain.
  if (!isStronglyConnected(P)) {
    return;
  }
  const double sqrtEps = std::sqrt(std::numeric_limits<double>::epsilon());

  // The eigenvalue computation is shared by slem, spectralGap and relaxationTime, and the
  // mixing time needs its own loop: do each of them only if requested.
  bool wantSpectral = false, wantMixing = false;
  for (int s = 0; s < k; ++s) {
    if (which[s] <= 2) wantSpectral = true;
    if (which[s] == 3) wantMixing = true;
  }

  // SLEM: the largest modulus among the eigenvalues other than the one
  // closest to 1, which must be within sqrt(eps) * max(1, n) of 1.
  double slem = NA_REAL;
  if (wantSpectral) {
    if (n == 1) {
      slem = 0;
    } else {
      arma::cx_vec values;
      if (arma::eig_gen(values, P)) {
        arma::uword unit = arma::index_min(arma::abs(values - 1.0));
        double tol = sqrtEps * std::max(1, n);
        if (std::abs(values(unit) - 1.0) <= tol) {
          double best = 0;
          for (arma::uword i = 0; i < values.n_elem; ++i) {
            if (i != unit) best = std::max(best, std::abs(values(i)));
          }
          if (std::isfinite(best)) slem = best;
        }
      }
    }
  }

  // Mixing time: the first t with max_i TV(P^t[i, ], pi) <= epsilon, for an
  // aperiodic chain; NA beyond maxIter.
  double mixing = NA_REAL;
  if (wantMixing) {
    NumericMatrix Pr = wrap(P);
    if (periodKernel(Pr) == 1) {
      if (n == 1) {
        mixing = 0;
      } else {
        arma::vec pi = steadyStateErgodicMatrix(P);
        bool valid = pi.is_finite() && arma::accu(pi) > 0 && pi.min() >= -sqrtEps;
        if (valid) {
          pi.elem(arma::find(pi < 0)).zeros();
          pi /= arma::accu(pi);
          arma::mat Pt = P;                      // P^t
          for (int t = 1; ; ++t) {
            // worst = max over the starting states of the total variation distance
            // between row i of P^t and pi
            double worst = 0;
            for (int i = 0; i < n; ++i) {
              double d = 0;
              for (int j = 0; j < n; ++j) d += std::abs(Pt(i, j) - pi(j));
              worst = std::max(worst, 0.5 * d);
            }
            if (worst <= epsilon) {
              mixing = t;
              break;
            }
            if (t >= maxIter) break;
            Pt = Pt * P;
          }
        }
      }
    }
  }

  // Derive the requested statistics: spectral gap = 1 - slem, relaxation time = 1 / gap
  // (infinite if the gap is numerically zero).
  for (int s = 0; s < k; ++s) {
    switch (which[s]) {
      case 0: out[s] = slem; break;
      case 1: out[s] = ISNA(slem) ? NA_REAL : 1 - slem; break;
      case 2:
        if (!ISNA(slem)) {
          double gap = 1 - slem;
          out[s] = gap <= sqrtEps ? R_PosInf : 1 / gap;
        }
        break;
      case 3: out[s] = mixing; break;
    }
  }
}

// One bootstrap replicate of the transition counts, written into N.
// Multinomial: for each state i in turn, departures[i] draws from row i
// (normalized as stats::rmultinom does). Parametric: one path of each
// observed length from each observed first state.
void drawCounts(const arma::mat& P, const arma::mat& cum, const IntegerVector& departures,
                bool parametric, const IntegerVector& lengths, const IntegerVector& starts,
                arma::mat& N, std::vector<double>& prob, std::vector<int>& draw) {
  const int n = P.n_rows;
  N.zeros();
  if (!parametric) {
    // multinomial: row i of the new counts is Multinomial(departures[i], row i of P)
    for (int i = 0; i < n; ++i) {
      double sum = 0;
      for (int j = 0; j < n; ++j) if (P(i, j) > 0) sum += P(i, j);
      for (int j = 0; j < n; ++j) prob[j] = P(i, j) / sum;
      R::rmultinom(departures[i], prob.data(), n, draw.data());
      for (int j = 0; j < n; ++j) N(i, j) = draw[j];
    }
  } else {
    // parametric: simulate again every observed sequence (same first state and length) and
    // count its transitions
    for (int s = 0; s < lengths.size(); ++s) {
      int x = starts[s];
      for (int t = 1; t < lengths[s]; ++t) {
        int y = drawFromCumulative(cum.colptr(0) + x, n, n);
        N(x, y) += 1;
        x = y;
      }
    }
  }
}

// Validate the arguments shared by the exported functions.
arma::mat checkedMatrix(const NumericMatrix& Pr) {
  if (Pr.nrow() != Pr.ncol() || Pr.nrow() < 1) {
    stop("The transition matrix must be square.");
  }
  arma::mat P(const_cast<double*>(Pr.begin()), Pr.nrow(), Pr.ncol(), true);
  if (!P.is_finite() || P.min() < 0) {
    stop("The transition matrix must be finite and non-negative.");
  }
  return P;
}

} // namespace

// Statistics of one row-stochastic matrix (used by the tests to compare
// with the R functions).
// [[Rcpp::export(.spectralStatisticsRcpp)]]
NumericVector spectralStatisticsRcpp(NumericMatrix Pr, IntegerVector which,
                                     double epsilon, int maxIter) {
  arma::mat P = checkedMatrix(Pr);
  NumericVector out(which.size());
  spectralStatistics(P, which, epsilon, maxIter, out.begin());
  return out;
}

// One replicate of transition counts (the R path of bootstrapProperties(),
// used when some statistic is a user function).
// [[Rcpp::export(.bootstrapCountsRcpp)]]
NumericMatrix bootstrapCountsRcpp(NumericMatrix Pr, IntegerVector departures, bool parametric,
                                  IntegerVector lengths, IntegerVector starts) {
  arma::mat P = checkedMatrix(Pr);
  const int n = P.n_rows;
  arma::mat cum = rowCumulative(P);
  arma::mat N(n, n);
  std::vector<double> prob(n);
  std::vector<int> draw(n);
  drawCounts(P, cum, departures, parametric, lengths, starts, N, prob, draw);
  NumericMatrix out = wrap(N);
  out.attr("dimnames") = Pr.attr("dimnames");
  return out;
}

// All the replicates of the built-in statistics: one row per replicate.
// A parametric replicate in which some state is never left gives a row of
// NA, as in the R code.
// [[Rcpp::export(.bootstrapSpectralRcpp)]]
NumericMatrix bootstrapSpectralRcpp(NumericMatrix Pr, IntegerVector departures, int nboot,
                                    bool parametric, IntegerVector lengths, IntegerVector starts,
                                    IntegerVector which, double epsilon, int maxIter) {
  arma::mat P = checkedMatrix(Pr);
  const int n = P.n_rows;
  const int k = which.size();
  arma::mat cum = rowCumulative(P);
  arma::mat N(n, n), Pstar(n, n);
  std::vector<double> prob(n), row(k);
  std::vector<int> draw(n);
  NumericMatrix out(nboot, k);
  for (int b = 0; b < nboot; ++b) {
    if (b % 64 == 0) Rcpp::checkUserInterrupt();
    // new counts, then the chain estimated from them (row-normalized counts)
    drawCounts(P, cum, departures, parametric, lengths, starts, N, prob, draw);
    arma::vec totals = arma::sum(N, 1);
    if (totals.min() <= 0) {
      for (int s = 0; s < k; ++s) out(b, s) = NA_REAL;
      continue;
    }
    Pstar = N.each_col() / totals;
    spectralStatistics(Pstar, which, epsilon, maxIter, row.data());
    for (int s = 0; s < k; ++s) out(b, s) = row[s];
  }
  return out;
}

// Simulated sequences for plotSequence(): one row per sequence, the states
// as 1-based codes. The initial state is drawn from mu, the following ones
// from the rows of P, by inversion with R's uniform generator.
// [[Rcpp::export(.simulateCodesRcpp)]]
IntegerMatrix simulateCodesRcpp(NumericMatrix Pr, NumericVector mu, int nSequences, int steps) {
  arma::mat P = checkedMatrix(Pr);
  const int n = P.n_rows;
  if (mu.size() != n) {
    stop("mu must have one probability per state.");
  }
  // cum: cumulative probabilities of each row, for sampling by inversion of one uniform draw
  arma::mat cum = rowCumulative(P);
  std::vector<double> cumMu(n);
  double acc = 0;
  for (int j = 0; j < n; ++j) {
    acc += mu[j];
    cumMu[j] = acc;
  }
  IntegerMatrix out(nSequences, steps);
  for (int s = 0; s < nSequences; ++s) {
    if (s % 256 == 0) Rcpp::checkUserInterrupt();
    int x = drawFromCumulative(cumMu.data(), n, 1);
    out(s, 0) = x + 1;
    for (int t = 1; t < steps; ++t) {
      x = drawFromCumulative(cum.colptr(0) + x, n, n);
      out(s, t) = x + 1;
    }
  }
  return out;
}
