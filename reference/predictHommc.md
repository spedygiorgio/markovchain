# Simulate from a higher order multivariate Markov chain

Simulates the next `t` states of each sequence of a higher order
multivariate Markov chain, given their last `order` states.

## Usage

``` r
predictHommc(hommc, t, init)
```

## Arguments

- hommc:

  An [`hommc`](hommc-class.md) object, for instance fitted by
  [`fitHighOrderMultivarMC()`](fitHighOrderMultivarMC.md).

- t:

  Number of steps to simulate, a positive integer.

- init:

  A matrix with one row per sequence and `order` columns, the last
  states of each sequence, the most recent in the last column. If
  missing, every sequence starts from the first state, repeated.

## Value

A character matrix with one row per sequence and `t` columns, the
simulated states.

## Details

At each step the next state of sequence \\j\\ is drawn from \\\sum_k
\sum_h \lambda\_{jkh} P_h^{(jk)} x^{(k)}\_{t-h+1}\\, where
\\x^{(k)}\_{t-h+1}\\ is the state of sequence \\k\\ \\h - 1\\ steps
before the current one (Ching et al., 2008). The matrices are read by
column (`P[to, from]`), as returned by
[`fitHighOrderMultivarMC()`](fitHighOrderMultivarMC.md), or by row when
the `byrow` slot of `hommc` is `TRUE`; all the sequences are drawn from
the same past before it is updated.

## References

Ching, W.-K., Ng, M. K. and Fung, E. S. (2008). Higher-order
multivariate Markov chains and their applications. Linear Algebra and
its Applications 428(2-3), 492-507.
[doi:10.1016/j.laa.2007.05.021](https://doi.org/10.1016/j.laa.2007.05.021)

## See also

Other higher order chains:
[`fitHighOrderMultivarMC()`](fitHighOrderMultivarMC.md),
[`fitHigherOrder()`](fitHigherOrder.md)

## Author

Vandit Jain

## Examples

``` r
data <- matrix(c("2", "1", "3", "3", "4", "3", "2", "1", "3", "3", "2", "1",
                 "2", "4", "4", "4", "4", "2", "3", "3", "1", "4", "3", "3"),
               ncol = 2)
if (requireNamespace("Rsolnp", quietly = TRUE)) {
  fit <- fitHighOrderMultivarMC(data, order = 2, Norm = 2)
  # the last two states of each sequence
  init <- t(data[11:12, ])
  predictHommc(fit, 5, init)
}
#> This function is experimental
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,] "1"  "3"  "3"  "4"  "2" 
#> [2,] "1"  "4"  "4"  "2"  "3" 
```
