# Second largest eigenvalue modulus (SLEM) of a Markov chain

Computes the second largest eigenvalue modulus (SLEM) of a finite,
irreducible discrete-time Markov chain.

## Usage

``` r
slem(object)

# S4 method for class 'markovchain'
slem(object)
```

## Arguments

- object:

  A `markovchain` object representing a finite, irreducible
  discrete-time Markov chain.

## Value

A numeric scalar in \\\[0,1\]\\ containing the SLEM. For the trivial
one-state chain, `0` is returned.

## Details

For a row-stochastic transition matrix \\P\\, let
\\1=\lambda_1,\lambda_2,\ldots,\lambda_n\\ be its eigenvalues. By the
Perron-Frobenius theorem an irreducible chain has \\\lambda_1=1\\ with
algebraic multiplicity one, and \\\|\lambda_k\|\le 1\\ for every \\k\\.
The SLEM is \$\$\mathrm{SLEM} = \max\_{k\>1} \|\lambda_k\|.\$\$

Only irreducibility is required, not aperiodicity. If the chain is
periodic, at least one non-trivial eigenvalue also has modulus one (e.g.
\\\lambda=-1\\ for a 2-cycle), so `slem()` correctly returns `1` rather
than rejecting the chain: a periodic chain genuinely does not contract
towards its stationary distribution, which `SLEM = 1` reflects.

Repeated or complex non-trivial eigenvalues are handled through their
modulus [`Mod()`](https://rdrr.io/r/base/complex.html), so
complex-conjugate pairs contribute the same value and ties do not need
to be broken.

The implementation calls [`eigen()`](https://rdrr.io/r/base/eigen.html)
with `only.values = TRUE`, so it never computes eigenvectors. Its time
complexity is \\O(n^3)\\ and its memory use is \\O(n^2)\\ for a dense
\\n\\-state transition matrix. It supports both row- and
column-stochastic storage.

## References

Levin, D. A. and Peres, Y. (2017). *Markov Chains and Mixing Times*, 2nd
edition. American Mathematical Society.

## See also

[`spectralGap`](spectralGap.md),
[`impliedTimescales`](impliedTimescales.md),
[`is.irreducible`](is.irreducible.md), [`period`](structuralAnalysis.md)

## Examples

``` r
statesNames <- c("a", "b")
mc <- new("markovchain",
  states = statesNames,
  transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
    byrow = TRUE, nrow = 2,
    dimnames = list(statesNames, statesNames)))
slem(mc)
#> [1] 0.6
```
