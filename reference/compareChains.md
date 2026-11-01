# Compare two Markov chains

Computes distances and divergences between two Markov chains on the same
states and, when the chains are given as observed sequences, tests
whether the sequences come from the same chain.

## Usage

``` r
compareChains(
  x,
  y,
  base = exp(1),
  test = c("G", "Pearson", "simulation", "none"),
  ...
)
```

## Arguments

- x, y:

  Two `markovchain` objects with the same states (matched by name), or
  two observed sequences (character, factor or numeric vectors), from
  which the chains are estimated by maximum likelihood on the union of
  their states.

- base:

  The base of the logarithm used by the divergences (default `exp(1)`,
  nats).

- test:

  For sequences, the test of equality of the two transition matrices
  passed to [`verifyHomogeneity`](statisticalTests.md): `"G"`
  (likelihood ratio, the default), `"Pearson"`, `"simulation"` (Monte
  Carlo), or `"none"`. Ignored for `markovchain` objects, whose sampling
  variability is not known.

- ...:

  Further arguments for [`verifyHomogeneity`](statisticalTests.md), such
  as `B` and `seed`.

## Value

An object of class `"markovchainComparison"`, a list with:

- `distances`:

  A named vector: `totalVariation`, the largest row-wise total variation
  distance \\\max_i \frac12\sum_j \|P\_{ij} - Q\_{ij}\|\\;
  `weightedTotalVariation`, the same row distances averaged with the
  stationary distribution of `x`; `hellinger`, the largest row-wise
  Hellinger distance \\\max_i \sqrt{\frac12 \sum_j (\sqrt{P\_{ij}} -
  \sqrt{Q\_{ij}})^2}\\; `frobenius`, \\\\P - Q\\\_F\\; `maximum`,
  \\\max\_{i,j}\|P\_{ij} - Q\_{ij}\|\\; and `stationary`, the total
  variation distance between the two stationary distributions.

- `divergences`:

  A named vector: `rateXY` and `rateYX`, the Kullback-Leibler divergence
  rates \\D(P\\Q)\\ and \\D(Q\\P)\\ (see
  [`divergenceRate`](divergenceRate.md)); `maxRowXY`, the largest
  row-wise Kullback-Leibler divergence of \\P\\ from \\Q\\; and
  `jensenShannon`, the largest row-wise Jensen-Shannon divergence, which
  is symmetric and finite.

- `stationaryDistributions`:

  A two-row matrix with the stationary distributions of `x` and `y`.

- `test`:

  For sequences, the `htest` returned by
  [`verifyHomogeneity`](statisticalTests.md); otherwise `NULL`.

- `x`, `y`:

  The two chains, as `markovchain` objects.

## Details

The measures that depend on a stationary distribution
(`weightedTotalVariation`, `stationary`, `rateXY`, `rateYX`) are `NA`
when the chain concerned has more than one. The row-wise measures and
their names follow PyDTMC's `calculate_transition_distance`,
`calculate_transition_divergence` and `calculate_stationary_distance`,
which give the same values.

For two sequences the test is that of Anderson and Goodman (1957) and
Kullback, Kupperman and Ku (1962): the hypothesis that both come from
the same chain, against different transition matrices. Each sequence is
treated as one sample; for several sequences per group, or for counts,
call [`verifyHomogeneity`](statisticalTests.md) directly.

## References

Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about
Markov chains. *The Annals of Mathematical Statistics*, 28(1), 89-110.

Kullback, S., Kupperman, M. and Ku, H. H. (1962). Tests for contingency
tables and Markov chains. *Technometrics*, 4(4), 573-608.

Rached, Z., Alajaji, F. and Campbell, L. L. (2004). The Kullback-Leibler
divergence rate between Markov sources. *IEEE Transactions on
Information Theory*, 50(5), 917-921.
[doi:10.1109/TIT.2004.826687](https://doi.org/10.1109/TIT.2004.826687)

## See also

[`divergenceRate`](divergenceRate.md),
[`verifyHomogeneity`](statisticalTests.md),
[`autoplot.markovchain`](autoplot.markovchain.md)
(`type = "comparison"`)

## Examples

``` r
s <- c("a", "b", "c")
P <- new("markovchain", states = s, transitionMatrix = matrix(
  c(0.5, 0.3, 0.2, 0.1, 0.6, 0.3, 0.25, 0.25, 0.5), 3, byrow = TRUE,
  dimnames = list(s, s)))
Q <- new("markovchain", states = s, transitionMatrix = matrix(
  c(0.4, 0.4, 0.2, 0.2, 0.5, 0.3, 0.3, 0.1, 0.6), 3, byrow = TRUE,
  dimnames = list(s, s)))
compareChains(P, Q)
#> Comparison of two Markov chains on 3 states
#> 
#> Distances:
#>         totalVariation weightedTotalVariation              hellinger 
#>                0.15000                0.11720                0.14250 
#>              frobenius                maximum             stationary 
#>                0.27390                0.15000                0.08825 
#> 
#> Divergences (rates weighted by the stationary distributions):
#>        rateXY        rateYX      maxRowXY jensenShannon 
#>       0.05428       0.05066       0.09233       0.02002 

# two sequences: distances between the estimated chains, and a test
set.seed(1)
x <- rmarkovchain(300, P, t0 = "a")
y <- rmarkovchain(300, Q, t0 = "a")
compareChains(x, y)
#> Comparison of two Markov chains on 3 states
#> 
#> Distances:
#>         totalVariation weightedTotalVariation              hellinger 
#>                0.22240                0.16940                0.16930 
#>              frobenius                maximum             stationary 
#>                0.38190                0.22240                0.04333 
#> 
#> Divergences (rates weighted by the stationary distributions):
#>        rateXY        rateYX      maxRowXY jensenShannon 
#>       0.09151       0.09950       0.10150       0.02804 
#> 
#> Test that both sequences come from the same chain:
#>   G-squared = 27.89, df = 6, p-value = 9.851e-05
```
