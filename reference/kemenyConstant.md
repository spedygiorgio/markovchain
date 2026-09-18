# Kemeny's constant of a Markov chain

Computes Kemeny's constant for a finite, irreducible discrete-time
Markov chain. It is the stationary-distribution-weighted mean hitting
time of a randomly selected destination and is independent of the
starting state.

## Usage

``` r
kemenyConstant(object)

# S4 method for class 'markovchain'
kemenyConstant(object)
```

## Arguments

- object:

  A `markovchain` object representing a finite, irreducible
  discrete-time Markov chain.

## Value

A numeric scalar containing Kemeny's constant.

## Details

For a row-stochastic transition matrix \\P\\, let \\\pi\\ be its unique
stationary distribution and define \$\$Z = (I - P +
\mathbf{1}\pi^T)^{-1}.\$\$ With hitting times defined by \\T_j = \inf\\n
\ge 0: X_n=j\\\\, so that \\m\_{jj}=0\\, the function returns \$\$K =
\sum_j \pi_j m\_{ij} = \mathrm{tr}(Z)-1.\$\$ The value does not depend
on the starting state \\i\\.

Irreducibility is sufficient; aperiodicity is not required. Reducible
chains can have multiple stationary distributions and are rejected.

Some references instead put the mean first-return time
\\m\_{jj}=1/\pi_j\\ on the diagonal. Under that convention the
corresponding stationary weighted sum is \\K+1\\, not \\K\\. This
function uses the zero-diagonal hitting-time convention, consistently
with [`meanFirstPassageTime()`](meanFirstPassageTime.md).

The implementation uses a dense LAPACK solve for the fundamental matrix
\\Z\\. Its time complexity is \\O(n^3)\\ and its memory use is
\\O(n^2)\\, as expected for a dense exact computation. It supports both
row- and column-stochastic storage.

## References

Kemeny, J. G. and Snell, J. L. (1960). *Finite Markov Chains*. D. Van
Nostrand, Princeton, NJ.

## See also

[`meanFirstPassageTime`](meanFirstPassageTime.md),
[`steadyStates`](steadyStates.md), [`is.irreducible`](is.irreducible.md)

## Examples

``` r
statesNames <- c("a", "b")
mc <- new("markovchain",
  states = statesNames,
  transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
    byrow = TRUE, nrow = 2,
    dimnames = list(statesNames, statesNames)))
kemenyConstant(mc)
#> [1] 2.5
```
