# Fundamental matrix of an absorbing Markov chain

Computes the fundamental matrix of an absorbing discrete-time Markov
chain. For the transient-state submatrix \\Q\\, the fundamental matrix
is \\N = (I - Q)^{-1}\\. Its \\(i,j)\\ entry is the expected number of
visits to transient state \\j\\ before absorption when starting from
transient state \\i\\.

## Usage

``` r
fundamentalMatrix(object)
```

## Arguments

- object:

  A `markovchain` object representing an absorbing Markov chain.

## Value

A numeric matrix containing the fundamental matrix, with transient state
names as row and column names.

## Details

An absorbing Markov chain must contain at least one absorbing state, and
all recurrent states must be absorbing. The function extracts the
transition submatrix corresponding to transient states and computes its
inverse complement. If there are no transient states, an empty matrix is
returned.

## References

Kemeny, J. G. and Snell, J. L. (1976). \*Finite Markov Chains\*.
Springer.

## See also

[`absorbingStates`](structuralAnalysis.md),
[`transientStates`](structuralAnalysis.md),
[`meanAbsorptionTime`](meanAbsorptionTime.md),
[`absorptionProbabilities`](absorptionProbabilities.md)

## Examples

``` r
states <- c("a", "b", "absorbed")
mc <- new("markovchain", states = states,
  transitionMatrix = matrix(c(
    0.5, 0.4, 0.1,
    0.2, 0.6, 0.2,
    0,   0,   1
  ), nrow = 3, byrow = TRUE,
  dimnames = list(states, states)))

fundamentalMatrix(mc)
#>          a        b
#> a 3.333333 3.333333
#> b 1.666667 4.166667
```
