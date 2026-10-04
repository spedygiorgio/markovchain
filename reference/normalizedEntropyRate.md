# Normalized entropy rate of a Markov chain

The entropy rate of a finite irreducible Markov chain divided by the
topological entropy of its graph, a measure of how random the chain is
relative to the most random chain with the same possible transitions.

## Usage

``` r
normalizedEntropyRate(object)

# S4 method for class 'markovchain'
normalizedEntropyRate(object)
```

## Arguments

- object:

  A `markovchain` object representing a finite, irreducible
  discrete-time Markov chain.

## Value

A numeric scalar in \\\[0,1\]\\: `0` for deterministic dynamics, `1`
when the entropy rate reaches the topological entropy.

## Details

The ratio \\H / h\_{top}\\ does not depend on the logarithm base, so no
`base` argument is needed. By the variational principle
([`topologicalEntropy`](topologicalEntropy.md)) it never exceeds one;
tiny excursions due to round-off are clipped to \\\[0,1\]\\. When
\\h\_{top}=0\\ the chain has a single possible path, hence entropy rate
zero, and the ratio is defined here to be `0`, the convention used by
PyDTMC's `entropy_rate_normalized`.

Irreducibility is required, as in [`entropyRate`](entropyRate.md), which
is what makes the entropy rate well defined.

## See also

[`entropyRate`](entropyRate.md),
[`topologicalEntropy`](topologicalEntropy.md)

## Examples

``` r
statesNames <- c("a", "b")
mc <- new("markovchain",
  states = statesNames,
  transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
    byrow = TRUE, nrow = 2,
    dimnames = list(statesNames, statesNames)))
normalizedEntropyRate(mc)
#> [1] 0.5720694
```
