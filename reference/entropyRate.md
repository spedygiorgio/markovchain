# Entropy rate of a Markov chain

Computes the entropy rate of a finite, irreducible discrete-time Markov
chain from its stationary distribution and transition matrix.

## Usage

``` r
entropyRate(object, base = 2)

# S4 method for class 'markovchain'
entropyRate(object, base = 2)
```

## Arguments

- object:

  A `markovchain` object representing a finite, irreducible
  discrete-time Markov chain.

- base:

  A finite numeric scalar strictly greater than one. The default, `2`,
  returns entropy in bits per transition. Use `exp(1)` for nats per
  transition.

## Value

A non-negative numeric scalar containing the entropy rate in units
determined by `base`.

## Details

For a row-stochastic transition matrix \\P\\ and stationary distribution
\\\pi\\, the entropy rate is \$\$H = -\sum_i \pi_i \sum_j
p\_{ij}\log_b(p\_{ij}),\$\$ with zero-probability transitions
contributing zero by continuity.

For a stationary first-order Markov chain, the entropy rate equals the
conditional entropy \\H(X\_{t+1}\mid X_t)\\. Irreducibility guarantees a
unique stationary distribution; aperiodicity is not required.

Reducible chains may admit multiple stationary distributions and
therefore different entropy rates. This method rejects them rather than
silently selecting one stationary distribution.

Transitions with probability zero are ignored, implementing the standard
convention \\0\log(0)=0\\ without evaluating `log(0)`.

The stationary-distribution computation dominates the running time. Once
the stationary distribution is available, evaluating the entropy rate
takes \\O(n^2)\\ time and \\O(n^2)\\ temporary memory for a dense
\\n\\-state transition matrix.

## References

Cover, T. M. and Thomas, J. A. (2006). *Elements of Information Theory*,
2nd edition. Wiley.

Strelioff, C. C., Crutchfield, J. P. and Huebler, A. W. (2007).
Inferring Markov chains: Bayesian estimation, model comparison, entropy
rate, and out-of-class modeling. *Physical Review E*, 76, 011106.

## See also

[`steadyStates`](steadyStates.md), [`is.irreducible`](is.irreducible.md)

## Examples

``` r
statesNames <- c("a", "b")
mc <- new("markovchain",
  states = statesNames,
  transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
    byrow = TRUE, nrow = 2,
    dimnames = list(statesNames, statesNames)))
entropyRate(mc)
#> [1] 0.5720694
entropyRate(mc, base = exp(1))
#> [1] 0.3965283
```
