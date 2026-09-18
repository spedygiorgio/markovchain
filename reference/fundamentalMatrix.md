# Fundamental matrix of an absorbing Markov chain

Computes the fundamental matrix of a finite absorbing discrete-time
Markov chain. If \\Q\\ is the transition submatrix restricted to
transient states, the fundamental matrix is \$\$N = I + Q + Q^2 + \cdots
= (I - Q)^{-1}.\$\$ The \\(i,j)\\ entry is the expected number of visits
to transient state \\j\\, including the initial visit when \\i = j\\,
before absorption, when the chain starts in transient state \\i\\.

## Usage

``` r
fundamentalMatrix(object)
```

## Arguments

- object:

  A `markovchain` object representing an absorbing Markov chain.

## Value

A numeric matrix containing the fundamental matrix, with transient state
names as row and column names. If all states are absorbing, the result
is a 0-by-0 matrix because there are no transient states.

## Details

For a finite absorbing chain, the state space can be reordered so that
the transition matrix has canonical form with transient block \\Q\\ and
an absorbing block. The spectral radius of \\Q\\ is less than one, so
the Neumann series converges and \\I-Q\\ is nonsingular.

The fundamental matrix also gives the expected time to absorption
through \\t = N 1\\, where \\1\\ is a vector of ones, and absorption
probabilities through \\B = N R\\, where \\R\\ contains transition
probabilities from transient to absorbing states.

The function requires an absorbing Markov chain: at least one absorbing
state must exist and every recurrent state must be absorbing. A chain
with no transient states is a valid degenerate case and returns a 0-by-0
matrix.

## References

Kemeny, J. G. and Snell, J. L. (1976). \*Finite Markov Chains\*.
Springer.

Grinstead, C. M. and Snell, J. L. (1997). \*Introduction to
Probability\*. American Mathematical Society.

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
