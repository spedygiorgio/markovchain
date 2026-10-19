// mtdFreeWeights.cpp: maximum likelihood of the MTD model with lag weights of
// any sign, the core of fitMTD(negativeWeights = TRUE) (R/fitMTD.R,
// .mtdFreeWeights()).
//
// For a junior reader. The model gives the probability of moving to state j
// after the states i_1, ..., i_k as f = sum_g lambda_g Q[i_g, j]. With
// negative weights some f could become negative, so the parameters must
// satisfy Raftery's (1985) condition. The likelihood is maximized by block
// coordinate ascent: in each block that sums to one (the weights, then each
// row of Q seen in the data) a little mass moves from the entry with the
// smallest partial derivative to the one with the largest, by the step that
// maximizes the likelihood along that direction within the feasible
// interval. Both the largest feasible step and the best step are found by
// bisection (60 halvings each). In R each move took hundreds of small
// function calls, so a fit could take minutes; the algorithm here is the same,
// step by step, with the sums accumulated in long double as R's sum() and
// colSums() do. Results agree with the R code up to rounding (the tests
// compare them).
//
// Data: `patterns` has one row per distinct observed pattern, the arrival
// state in column 0 and the departure state at lag g in column g (0-based
// codes); `counts` their numbers of occurrences.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// Floating-point contraction off: clang (the compiler of CRAN's and GitHub's
// macOS arm64 builds) by default fuses a * b + c into one fused multiply-add,
// rounded once, while R rounds the product and the sum separately. This keeps
// the operations those of the R code; the results are bit for bit R's on
// x86-64, and may still differ in the last bits elsewhere (they do on macOS
// arm64, where long double is double).
#ifdef __clang__
#pragma clang fp contract(off)
#endif

using namespace Rcpp;

namespace {

struct MtdProblem {
  int P;                        // number of patterns
  int k;                        // number of lags
  int r;                        // number of states
  std::vector<int> to;          // arrival state of each pattern
  std::vector<int> from;        // from[p + P * g]: departure state at lag g
  std::vector<double> counts;
  std::vector<int> rows;        // rows of Q seen in the data (0-based)
};

// Raftery's condition, as .mtdFeasible(): all of Q non-negative, and for
// every arrival state j, min_i Q[i, j] * (sum of positive weights) -
// max_i Q[i, j] * (sum of absolute negative weights) >= -tol, i in rows.
bool feasible(const MtdProblem& pb, const std::vector<double>& lambda,
              const std::vector<double>& Q, double tol = 1e-12) {
  const int r = pb.r;
  for (double q : Q) {
    if (!(q >= 0)) return false;
  }
  // sum(pmax(lambda, 0)) and sum(pmax(-lambda, 0)), accumulated in long
  // double as R's sum() does: after a move that stops on the boundary of the
  // condition the slack sits on the tolerance, and the last bits decide.
  long double posL = 0, negL = 0;
  for (double l : lambda) {
    posL += l > 0 ? l : 0;
    negL += -l > 0 ? -l : 0;
  }
  const double pos = (double) posL, neg = (double) negL;
  double slack = std::numeric_limits<double>::infinity();
  for (int j = 0; j < r; ++j) {
    double lo = std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();
    for (int i : pb.rows) {
      double q = Q[i + r * j];
      if (q < lo) lo = q;
      if (q > hi) hi = q;
    }
    double s = lo * pos - hi * neg;
    if (s < slack) slack = s;
  }
  return slack >= -tol;
}

// f_p = sum_g lambda_g Q[from(p, g), to(p)]
void mixture(const MtdProblem& pb, const std::vector<double>& lambda,
             const std::vector<double>& Q, std::vector<double>& f) {
  for (int p = 0; p < pb.P; ++p) {
    double s = 0;
    for (int g = 0; g < pb.k; ++g) {
      s += Q[pb.from[p + pb.P * g] + pb.r * pb.to[p]] * lambda[g];
    }
    f[p] = s;
  }
}

double logLik(const MtdProblem& pb, const std::vector<double>& f) {
  long double s = 0;
  for (int p = 0; p < pb.P; ++p) {
    if (!(f[p] > 0)) return R_NegInf;
    s += pb.counts[p] * std::log(f[p]);
  }
  return (double) s;
}

// f + delta h stays positive for delta below this bound
double positiveBound(const std::vector<double>& f, const std::vector<double>& h) {
  double best = R_PosInf;
  bool any = false;
  for (size_t p = 0; p < h.size(); ++p) {
    if (h[p] < 0) {
      any = true;
      double v = f[p] / -h[p];
      if (v < best) best = v;
    }
  }
  return any ? best * (1 - 1e-12) : R_PosInf;
}

// Largest feasible step in [0, upper] by bisection (the feasible set is
// convex), as maxStep() in R.
template <typename Feasible>
double maxStep(double upper, Feasible ok) {
  if (ok(upper)) return upper;
  double lo = 0, hi = upper;
  for (int it = 0; it < 60; ++it) {
    double mid = (lo + hi) / 2;
    if (ok(mid)) lo = mid; else hi = mid;
  }
  return lo;
}

// Maximizer of sum c log(f + delta h) on [0, upper] by bisection on the
// derivative (concave in delta), as lineMax() in R.
double lineMax(const MtdProblem& pb, const std::vector<double>& f,
               const std::vector<double>& h, double upper) {
  auto d1 = [&](double delta) {
    long double s = 0;
    for (int p = 0; p < pb.P; ++p) s += pb.counts[p] * h[p] / (f[p] + delta * h[p]);
    return (double) s;
  };
  if (d1(upper) >= 0) return upper;
  double lo = 0, hi = upper;
  for (int it = 0; it < 60; ++it) {
    double mid = (lo + hi) / 2;
    if (d1(mid) > 0) lo = mid; else hi = mid;
  }
  return lo;
}

// index of the first maximum / minimum, as which.max() / which.min()
int firstMax(const std::vector<double>& v) {
  int a = 0;
  for (size_t i = 1; i < v.size(); ++i) if (v[i] > v[a]) a = (int) i;
  return a;
}

} // namespace

// [[Rcpp::export(.mtdFreeWeightsRcpp)]]
List mtdFreeWeightsRcpp(IntegerMatrix patterns, NumericVector counts, NumericVector lambda0,
                        NumericMatrix Q0, IntegerVector rows, double tol, int maxit) {
  MtdProblem pb;
  pb.P = patterns.nrow();
  pb.k = patterns.ncol() - 1;
  pb.r = Q0.nrow();
  if (pb.k < 1 || (int) lambda0.size() != pb.k || Q0.ncol() != pb.r ||
      (int) counts.size() != pb.P || rows.size() < 1) {
    stop("inconsistent arguments.");
  }
  pb.to.resize(pb.P);
  pb.from.resize((size_t) pb.P * pb.k);
  pb.counts.assign(counts.begin(), counts.end());
  for (int p = 0; p < pb.P; ++p) {
    pb.to[p] = patterns(p, 0);
    for (int g = 0; g < pb.k; ++g) pb.from[p + pb.P * g] = patterns(p, g + 1);
  }
  for (int p = 0; p < pb.P * (pb.k + 1); ++p) {
    int code = patterns[p];
    if (code < 0 || code >= pb.r) stop("pattern codes must lie between 0 and r - 1.");
  }
  for (int i : rows) {
    if (i < 1 || i > pb.r) stop("rows must lie between 1 and r.");
    pb.rows.push_back(i - 1);
  }
  const int P = pb.P, k = pb.k, r = pb.r;
  std::vector<double> lambda(lambda0.begin(), lambda0.end());
  std::vector<double> Q(Q0.begin(), Q0.end());   // column-major, Q[i + r j]
  std::vector<double> f(P), h(P), w(P), lamI(P), grad, trial;

  mixture(pb, lambda, Q, f);
  double L = logLik(pb, f);
  if (!R_FINITE(L) || !feasible(pb, lambda, Q)) {
    stop("the starting point of the MTD fit with negative weights is not feasible");
  }

  int sweeps = 0;
  bool converged = false;
  while (sweeps < maxit) {
    ++sweeps;
    if (sweeps % 16 == 0) Rcpp::checkUserInterrupt();
    const double Lstart = L;

    // weights: grad_g = sum_p C[p, g] counts_p / f_p
    for (int p = 0; p < P; ++p) w[p] = pb.counts[p] / f[p];
    grad.assign(k, 0);
    for (int g = 0; g < k; ++g) {
      long double s = 0;
      for (int p = 0; p < P; ++p) s += Q[pb.from[p + P * g] + r * pb.to[p]] * w[p];
      grad[g] = (double) s;
    }
    {
      int a = firstMax(grad);
      int b = 0;
      for (int g = 1; g < k; ++g) if (grad[g] < grad[b]) b = g;
      if (grad[a] - grad[b] > 1e-12 * std::abs(L)) {
        for (int p = 0; p < P; ++p) {
          h[p] = Q[pb.from[p + P * a] + r * pb.to[p]] - Q[pb.from[p + P * b] + r * pb.to[p]];
        }
        double upper = std::min(positiveBound(f, h), 10.0);
        upper = maxStep(upper, [&](double d) {
          trial = lambda;
          trial[a] += d;
          trial[b] -= d;
          return feasible(pb, trial, Q);
        });
        if (upper > 0) {
          double delta = lineMax(pb, f, h, upper);
          lambda[a] += delta;
          lambda[b] -= delta;
          for (int p = 0; p < P; ++p) f[p] += delta * h[p];
        }
      }
    }

    // rows of Q
    for (int i : pb.rows) {
      for (int p = 0; p < P; ++p) {
        w[p] = pb.counts[p] / f[p];
        // weights of the lags at which the departure state is i
        double s = 0;
        for (int g = 0; g < k; ++g) if (pb.from[p + P * g] == i) s += lambda[g];
        lamI[p] = s;
      }
      grad.assign(r, 0);
      {
        std::vector<long double> acc(r, 0);
        for (int p = 0; p < P; ++p) acc[pb.to[p]] += w[p] * lamI[p];
        for (int j = 0; j < r; ++j) grad[j] = (double) acc[j];
      }
      int a = firstMax(grad);
      int b = -1;
      for (int j = 0; j < r; ++j) {
        if (Q[i + r * j] > 0 && (b < 0 || grad[j] < grad[b])) b = j;
      }
      if (b < 0 || grad[a] - grad[b] <= 1e-12 * std::abs(L)) continue;
      for (int p = 0; p < P; ++p) {
        h[p] = lamI[p] * ((pb.to[p] == a ? 1.0 : 0.0) - (pb.to[p] == b ? 1.0 : 0.0));
      }
      double upper = std::min(Q[i + r * b], positiveBound(f, h));
      upper = maxStep(upper, [&](double d) {
        trial = Q;
        trial[i + r * a] += d;
        trial[i + r * b] = std::max(trial[i + r * b] - d, 0.0);
        return feasible(pb, lambda, trial);
      });
      if (!(upper > 0)) continue;
      double delta = lineMax(pb, f, h, upper);
      Q[i + r * a] += delta;
      Q[i + r * b] = std::max(Q[i + r * b] - delta, 0.0);
      for (int p = 0; p < P; ++p) f[p] += delta * h[p];
    }

    mixture(pb, lambda, Q, f);
    L = logLik(pb, f);
    if (L - Lstart <= tol * (std::abs(L) + tol)) {
      converged = true;
      break;
    }
  }

  NumericMatrix Qout(r, r);
  std::copy(Q.begin(), Q.end(), Qout.begin());
  return List::create(_["lambda"] = NumericVector(lambda.begin(), lambda.end()),
                      _["Q"] = Qout, _["logLik"] = L, _["iterations"] = sweeps,
                      _["converged"] = converged);
}
