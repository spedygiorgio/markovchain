# Relaxation time of a Markov chain

The relaxation time \\1/\mathrm{gap}\\ of a finite irreducible
discrete-time Markov chain, the reciprocal of its spectral gap.

## Usage

``` r
relaxationTime(object)

# S4 method for class 'markovchain'
relaxationTime(object)
```

## Arguments

- object:

  A `markovchain` object representing a finite, irreducible
  discrete-time Markov chain.

## Value

A positive numeric scalar, in number of steps: `Inf` for a periodic
chain (spectral gap zero), `1` for the trivial one-state chain.

## Details

Following Levin and Peres (2017, Section 12.2) the relaxation time is
\\t\_{rel} = 1/\gamma\\, with \\\gamma = 1 - \mathrm{SLEM}\\ the
spectral gap of [`spectralGap`](spectralGap.md). It is the quantity
PyDTMC calls `relaxation_rate`, even if it is a time, not a rate. The
two are the same number.

PyDTMC's `mixing_rate` is \\-1/\log(\mathrm{SLEM})\\, i.e. the timescale
of the slowest non-trivial mode, which is already returned as the first
element (`"tau2"`) of [`impliedTimescales`](impliedTimescales.md). It is
therefore not repeated as a separate function.

## References

Levin, D. A. and Peres, Y. (2017). *Markov Chains and Mixing Times*, 2nd
edition. American Mathematical Society.

## See also

[`spectralGap`](spectralGap.md), [`slem`](slem.md),
[`impliedTimescales`](impliedTimescales.md),
[`mixingTime`](mixingTime.md)

## Examples

``` r
statesNames <- c("a", "b")
mc <- new("markovchain",
  states = statesNames,
  transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
    byrow = TRUE, nrow = 2,
    dimnames = list(statesNames, statesNames)))
relaxationTime(mc)
#> [1] 2.5
```
