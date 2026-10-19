// markovPropertyTest.cpp: the Monte Carlo part of verifyMarkovProperty(method =
// "simulation") (R/statisticalTests.R).
//
// For a junior reader. The test compares the counts of triples of
// consecutive states (i, j, k) with those expected if the next state depended
// only on the present one. With method = "simulation" the null distribution
// of the statistic is obtained by simulating B sequences from the estimated
// first-order chain and recomputing the statistic on each. In R each step of
// each sequence called sample.int(), which made the test take minutes with
// the default B = 9999.
//
// The functions here reproduce the R code exactly, not only in distribution:
// - one step draws the next state as sample.int(r, 1, prob = P[x, ]) does in
//   R (src/main/random.c: FixupProb, then ProbSampleNoReplace, which sorts the
//   probabilities with revsort() and inverts their cumulative sums), with R's
//   unif_rand(), so the same seed gives the same sequences;
// - the statistic is computed with the same operations in the same order,
//   the cells taken in the column-major order of the R array and summed in
//   long double as R's sum() does, so the same counts give the same value.
// The observed statistic is still computed in R; the tests compare the
// simulated statistics with those of the former R loop.

#include <Rcpp.h>
#include <R_ext/Utils.h>
#include <algorithm>
#include <cmath>
#include <vector>

// Floating-point contraction off: clang (the compiler of CRAN's and GitHub's
// macOS arm64 builds) by default fuses a * b + c into one fused multiply-add,
// rounded once, while R rounds the product and the sum separately; the last
// bits would then differ from R's. GCC on x86-64 never fuses without -mfma.
#ifdef __clang__
#pragma clang fp contract(off)
#endif

using namespace Rcpp;

namespace {

// Statistic of the Markov property test from the codes (0-based) of a
// sequence: counts N[i + j r + k r^2] of triples (past i, present j, future
// k), marginals N_ij and N_jk, expected counts N_ij N_jk / N_j, then the
// Pearson or G statistic over the cells with positive expectation (and, for
// G, positive count).
double markovStatistic(const int* x, int n, int r, bool pearson,
                       std::vector<int>& N, std::vector<int>& Nij,
                       std::vector<int>& Njk, std::vector<int>& Nj) {
  const int r2 = r * r;
  std::fill(N.begin(), N.end(), 0);
  for (int t = 0; t + 2 < n; ++t) {
    N[x[t] + x[t + 1] * r + x[t + 2] * r2] += 1;
  }
  std::fill(Nij.begin(), Nij.end(), 0);
  std::fill(Njk.begin(), Njk.end(), 0);
  std::fill(Nj.begin(), Nj.end(), 0);
  for (int k = 0; k < r; ++k) {
    for (int j = 0; j < r; ++j) {
      for (int i = 0; i < r; ++i) {
        int c = N[i + j * r + k * r2];
        Nij[i + j * r] += c;  // Nij[i, j]
        Njk[j + k * r] += c;  // Njk[j, k]
      }
    }
  }
  for (int k = 0; k < r; ++k) {
    for (int j = 0; j < r; ++j) Nj[j] += Njk[j + k * r];
  }
  long double sum = 0.0L;
  for (int k = 0; k < r; ++k) {
    for (int j = 0; j < r; ++j) {
      if (Nj[j] <= 0) continue;
      for (int i = 0; i < r; ++i) {
        // outer(Nij[, j], Njk[j, ]) / Nj[j]
        double e = (double) Nij[i + j * r] * (double) Njk[j + k * r] / (double) Nj[j];
        if (!(e > 0)) continue;
        double o = N[i + j * r + k * r2];
        if (pearson) {
          sum += (o - e) * (o - e) / e;
        } else if (o > 0) {
          sum += o * std::log(o / e);
        }
      }
    }
  }
  return pearson ? (double) sum : 2 * (double) sum;
}

} // namespace

// B statistics of sequences of length n simulated from the row-stochastic
// matrix P, each starting from the 1-based state x1.
// [[Rcpp::export(.markovPropertySimulationRcpp)]]
NumericVector markovPropertySimulationRcpp(NumericMatrix P, int x1, int n, int B, bool pearson) {
  const int r = P.nrow();
  if (P.ncol() != r || r < 1) stop("P must be a square matrix.");
  if (x1 < 1 || x1 > r) stop("x1 must lie between 1 and the number of states.");
  if (n < 3 || B < 1) stop("n must be at least 3 and B at least 1.");
  for (int i = 0; i < r; ++i) {
    double rowSum = 0;
    for (int j = 0; j < r; ++j) {
      double p = P(i, j);
      if (!R_FINITE(p) || p < 0) stop("P must be finite and non-negative.");
      rowSum += p;
    }
    if (!(rowSum > 0)) stop("every row of P must have a positive entry.");
  }
  std::vector<int> x(n), perm(r);
  std::vector<double> p(r);
  std::vector<int> N((size_t) r * r * r), Nij(r * r), Njk(r * r), Nj(r);
  NumericVector out(B);
  for (int b = 0; b < B; ++b) {
    if (b % 100 == 0) Rcpp::checkUserInterrupt();
    x[0] = x1 - 1;
    for (int t = 1; t < n; ++t) {
      const int i = x[t - 1];
      // FixupProb(): divide by the sum of the positive probabilities
      double sum = 0.0;
      for (int j = 0; j < r; ++j) {
        p[j] = P(i, j);
        if (p[j] > 0.0) sum += p[j];
      }
      for (int j = 0; j < r; ++j) {
        p[j] /= sum;
        perm[j] = j + 1;
      }
      // ProbSampleNoReplace() for one draw
      revsort(p.data(), perm.data(), r);
      double rT = unif_rand();
      double mass = 0;
      int j;
      for (j = 0; j < r - 1; ++j) {
        mass += p[j];
        if (rT <= mass) break;
      }
      x[t] = perm[j] - 1;
    }
    out[b] = markovStatistic(x.data(), n, r, pearson, N, Nij, Njk, Nj);
  }
  return out;
}
