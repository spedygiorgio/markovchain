# Generator matrix estimated from a matrix of transition frequencies

Estimates the generator of a continuous time Markov chain from the
empirical transition matrix observed over a period of length `t`, with
the regularization algorithms of Kreinin and Sidelnikova (2001), which
return a valid generator (non-negative off-diagonal rates, rows summing
to zero) even when the matrix logarithm of `P` is not one. The
computation is done by the suggested package ctmcd.

## Usage

``` r
freq2Generator(P, t = 1, method = "QO", logmethod = "Eigen")
```

## Arguments

- P:

  Matrix of relative transition frequencies over the period: a square
  matrix whose rows sum to one.

- t:

  Length of the observation period, a positive number (default 1).

- method:

  Regularization algorithm: `"QO"` (quasi-optimization, the default),
  `"WA"` (weighted adjustment) or `"DA"` (diagonal adjustment).

- logmethod:

  Method of the matrix logarithm used by ctmcd (default `"Eigen"`, see
  [`ctmcd::gmQO()`](https://rdrr.io/pkg/ctmcd/man/gmQO.html)).

## Value

The estimated generator, a matrix with the dimensions and dimnames of
`P`; `NULL`, with a warning, if ctmcd is not installed.

## References

Kreinin, A. and Sidelnikova, M. (2001). Regularization algorithms for
transition matrices. Algo Research Quarterly 4(1), 23-40.

## See also

Other continuous time Markov chains:
[`ExpectedTime()`](ExpectedTime.md), [`ctmcFit()`](ctmcFit.md),
[`generatorToTransitionMatrix()`](generatorToTransitionMatrix.md),
[`impreciseProbabilityatT()`](impreciseProbabilityatT.md),
[`is.CTMCirreducible()`](is.CTMCirreducible.md),
[`is.TimeReversible()`](is.TimeReversible.md),
[`probabilityatT()`](probabilityatT.md), [`rctmc()`](rctmc.md),
[`transition2Generator()`](transition2Generator.md)

## Examples

``` r
counts <- matrix(c(150, 2, 1, 1,
                   1, 200, 2, 1,
                   2, 1, 175, 1,
                   1, 1, 1, 150), nrow = 4, byrow = TRUE)
# relative frequencies, with the last state made absorbing
rel <- rbind((counts / rowSums(counts))[1:3, ], c(0, 0, 0, 1))
if (requireNamespace("ctmcd", quietly = TRUE)) freq2Generator(rel, 1)
#>              [,1]        [,2]         [,3] [,4]
#> [1,] -0.024212164  0.01544797  0.008764198    0
#> [2,]  0.006594821 -0.01822834  0.011633520    0
#> [3,]  0.013302567  0.00749703 -0.020799597    0
#> [4,]  0.000000000  0.00000000  0.000000000    0

# rating transitions (quasi-optimization estimate)
data(tm_abs)
tmRel <- rbind((tm_abs / rowSums(tm_abs))[1:7, ], c(rep(0, 7), 1))
if (requireNamespace("ctmcd", quietly = TRUE)) freq2Generator(tmRel, 1)
#>              AAA           AA            A          BBB           BB
#> AAA -0.109688198  0.104742772  0.004945426  0.000000000  0.000000000
#> AA   0.006375885 -0.095416840  0.088027218  0.001013738  0.000000000
#> A    0.000000000  0.037605369 -0.139127796  0.092863518  0.002082790
#> BBB  0.000000000  0.003101629  0.043766820 -0.100963213  0.044471253
#> BB   0.000000000  0.004024583  0.000000000  0.043976931 -0.142486437
#> B    0.000000000  0.005844572  0.003289642  0.005803757  0.058923108
#> C    0.000000000  0.000000000  0.000000000  0.000000000  0.006651241
#>      0.000000000  0.000000000  0.000000000  0.000000000  0.000000000
#>                 B            C           D
#> AAA  0.000000e+00  0.000000000 0.000000000
#> AA   0.000000e+00  0.000000000 0.000000000
#> A    1.064427e-05  0.004562575 0.002002900
#> BBB  4.257673e-03  0.001871779 0.003494059
#> BB   8.610403e-02  0.008380896 0.000000000
#> B   -1.932222e-01  0.064440303 0.054920850
#> C    1.547477e-01 -0.362361432 0.200962499
#>      0.000000e+00  0.000000000 0.000000000
```
