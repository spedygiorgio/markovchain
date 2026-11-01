# Generator matrix from a transition matrix

Computes the generator \\Q\\ of a continuous time Markov chain whose
transition matrix over a period of length `t` is `P`, as the matrix
logarithm \\Q = \log(P) / t\\, so that \\\exp(t Q) = P\\.

## Usage

``` r
transition2Generator(P, t = 1, method = "logarithm")
```

## Arguments

- P:

  A transition matrix (rows summing to one) over a period of length `t`.

- t:

  Length of the period, a positive number.

- method:

  Method used: only `"logarithm"`, the principal matrix logarithm
  computed by [`expm::logm()`](https://rdrr.io/pkg/expm/man/logm.html),
  is available.

## Value

A square matrix with the dimensions of `P`, the generator \\\log(P) /
t\\.

## Details

Not every transition matrix is the transition matrix of a continuous
time chain (the embedding problem): the matrix returned is a valid
generator only if its off-diagonal entries are non-negative, which is
not checked. [`freq2Generator()`](freq2Generator.md) gives regularized
estimates that are always valid generators.

## See also

Other continuous time Markov chains:
[`ExpectedTime()`](ExpectedTime.md), [`ctmcFit()`](ctmcFit.md),
[`freq2Generator()`](freq2Generator.md),
[`generatorToTransitionMatrix()`](generatorToTransitionMatrix.md),
[`impreciseProbabilityatT()`](impreciseProbabilityatT.md),
[`is.CTMCirreducible()`](is.CTMCirreducible.md),
[`is.TimeReversible()`](is.TimeReversible.md),
[`probabilityatT()`](probabilityatT.md), [`rctmc()`](rctmc.md)

## Examples

``` r
P <- matrix(c(0.4, 0.6, 0.1, 0.9), nrow = 2, byrow = TRUE)
Q <- transition2Generator(P)
Q
#>            [,1]       [,2]
#> [1,] -1.0319767  1.0319767
#> [2,]  0.1719961 -0.1719961
expm::expm(Q) # gives back P
#>      [,1] [,2]
#> [1,]  0.4  0.6
#> [2,]  0.1  0.9
```
