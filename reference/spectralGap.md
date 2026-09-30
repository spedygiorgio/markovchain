# Spectral gap of a Markov chain

Computes the spectral gap of a finite, irreducible discrete-time Markov
chain, a lightweight diagnostic of its convergence and mixing behaviour.

## Usage

``` r
spectralGap(object)

# S4 method for class 'markovchain'
spectralGap(object)
```

## Arguments

- object:

  A `markovchain` object representing a finite, irreducible
  discrete-time Markov chain.

## Value

A numeric scalar in \\\[0,1\]\\ containing the spectral gap. For the
trivial one-state chain, `1` is returned.

## Details

The spectral gap is defined from the second largest eigenvalue modulus
(SLEM, see [`slem`](slem.md)) as \$\$\mathrm{gap} = 1 -
\mathrm{SLEM}.\$\$

As with [`slem`](slem.md), only irreducibility is required. A periodic
chain has `SLEM = 1` and therefore spectral gap `0`: this is the
mathematically correct value, not an error condition, since a periodic
chain never contracts towards its stationary distribution. A larger
spectral gap indicates faster convergence to stationarity; see
[`impliedTimescales`](impliedTimescales.md) for the timescale associated
with each non-trivial eigenvalue individually, of which the SLEM gives
the slowest (dominant) one.

## References

Levin, D. A. and Peres, Y. (2017). *Markov Chains and Mixing Times*, 2nd
edition. American Mathematical Society.

## See also

[`slem`](slem.md), [`impliedTimescales`](impliedTimescales.md),
[`is.irreducible`](is.irreducible.md)

## Examples

``` r
statesNames <- c("a", "b")
mc <- new("markovchain",
  states = statesNames,
  transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
    byrow = TRUE, nrow = 2,
    dimnames = list(statesNames, statesNames)))
spectralGap(mc)
#> [1] 0.4
```
