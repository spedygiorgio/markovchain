# Kullback-Leibler divergence rate between two Markov chains

Computes the Kullback-Leibler divergence rate of a Markov chain `x` from
a Markov chain `y` on the same states: the limit of \\\frac{1}{n}
D(X_1,\ldots,X_n \\\\\\ Y_1,\ldots,Y_n)\\, the average information per
step for discriminating `x` against `y`.

## Usage

``` r
divergenceRate(x, y, weights = NULL, base = exp(1))
```

## Arguments

- x, y:

  Two `markovchain` objects with the same states, matched by name.

- weights:

  Optional non-negative numeric vector, one entry per state (named or in
  the order of `states(x)`), with which the row divergences are
  averaged; it is normalized to sum to one. The default, `NULL`, uses
  the stationary distribution of `x`, which must then be unique.

- base:

  The base of the logarithm: the default `exp(1)` gives the rate in
  nats, `2` in bits.

## Value

A single non-negative number, `Inf` when `x` makes with positive
probability a transition that `y` forbids from a state of positive
weight.

## Details

For chains with transition matrices \\P\\ (of `x`) and \\Q\\ (of `y`)
and with \\\pi\\ the stationary distribution of \\P\\, the divergence
rate is \$\$D(P \\\\\\ Q) = \sum_i \pi_i \sum_j P\_{ij} \log
\frac{P\_{ij}}{Q\_{ij}},\$\$ which Rached, Alajaji and Campbell (2004)
prove to be the limit above for irreducible chains, whatever the initial
distributions. By the ergodic theorem it is also the almost sure limit
of \\\frac1n \log \frac{\Pr_P(X_1, \ldots, X_n)}{\Pr_Q(X_1, \ldots,
X_n)}\\ along a path of `x`. It is zero exactly when the rows of \\P\\
and \\Q\\ agree on every state of positive weight, and it is not
symmetric.

With `weights` equal to the relative frequencies of the departure states
of a sequence, \\2N D(\hat P \\\\\\ Q)\\, with \\N\\ the number of
transitions, is the likelihood-ratio statistic \\G^2\\ of the hypothesis
that the sequence comes from \\Q\\ (Kullback, Kupperman and Ku 1962), as
computed by [`verifyEmpiricalToTheoretical`](statisticalTests.md).

## References

Rached, Z., Alajaji, F. and Campbell, L. L. (2004). The Kullback-Leibler
divergence rate between Markov sources. *IEEE Transactions on
Information Theory*, 50(5), 917-921.
[doi:10.1109/TIT.2004.826687](https://doi.org/10.1109/TIT.2004.826687)

Kullback, S., Kupperman, M. and Ku, H. H. (1962). Tests for contingency
tables and Markov chains. *Technometrics*, 4(4), 573-608.

## See also

[`compareChains`](compareChains.md), [`entropyRate`](entropyRate.md)

## Examples

``` r
s <- c("a", "b", "c")
P <- new("markovchain", states = s, transitionMatrix = matrix(
  c(0.5, 0.3, 0.2, 0.1, 0.6, 0.3, 0.25, 0.25, 0.5), 3, byrow = TRUE,
  dimnames = list(s, s)))
Q <- new("markovchain", states = s, transitionMatrix = matrix(
  c(0.4, 0.4, 0.2, 0.2, 0.5, 0.3, 0.3, 0.1, 0.6), 3, byrow = TRUE,
  dimnames = list(s, s)))
divergenceRate(P, Q)
#> [1] 0.05428363
divergenceRate(Q, P) # not symmetric
#> [1] 0.05065658
divergenceRate(P, P)
#> [1] 0
```
