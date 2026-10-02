# Functions to fit a higher order Markov chain

Given a sequence of states arising from a stationary state, it fits the
underlying Markov chain distribution with higher order.

## Usage

``` r
fitHigherOrder(sequence, order = 2, method = c("lsq", "mle"))
seq2freqProb(sequence)
seq2matHigh(sequence, order)
```

## Arguments

- sequence:

  A character list.

- order:

  Markov chain order

- method:

  How the weights \\\lambda\\ are estimated: `"lsq"` (default) or
  `"mle"`, see Details.

## Value

A list containing lambda, Q, and X.

## Details

The fitted model expresses the distribution of the next state as the
mixture \\\sum\_{i=1}^{k} \lambda_i Q_i x\_{t-i}\\ of the empirical
lag-\\i\\ transition matrices \\Q_i\\ (see `seq2matHigh`), with weights
\\\lambda_i \ge 0\\ summing to one. The matrices \\Q_i\\ are the same
for both methods; only the weights differ.

`method = "lsq"` (the default, and the only behaviour before the
argument existed) chooses \\\lambda\\ to minimize the squared distance
between the stationary distribution and its image under the mixture, as
in Ching et al.; it needs the Rsolnp package and returns `NULL` with a
message if it is unavailable.

`method = "mle"` chooses \\\lambda\\ to maximize the log-likelihood
\\\sum\_{t=k+1}^{n} \log \sum_i \lambda_i Q_i\[x_t, x\_{t-i}\]\\ of the
observations a model of order \\k\\ can predict. For fixed \\Q_i\\ the
problem is concave, so the maximum is global, and it is solved by the EM
algorithm for mixture weights, without Rsolnp. The weights are therefore
those that give the highest value of
[`higherOrderLogLik`](higherOrderLogLik.md) for the same observations.
Note that this is not the mixture transition distribution model of
Raftery (1985), in which a single matrix is shared by all lags and is
estimated together with the weights.

## References

Ching, W. K., Huang, X., Ng, M. K., & Siu, T. K. (2013). Higher-order
markov chains. In Markov Chains (pp. 141-176). Springer US.

Ching, W. K., Ng, M. K., & Fung, E. S. (2008). Higher-order multivariate
Markov chains and their applications. Linear Algebra and its
Applications, 428(2), 492-507.

Raftery, A. E. (1985). A model for high-order Markov chains. Journal of
the Royal Statistical Society, Series B, 47(3), 528-539.

## Author

Giorgio Spedicato, Tae Seung Kang

## Examples

``` r
sequence<-c("a", "a", "b", "b", "a", "c", "b", "a", "b", "c", "a", "b",
            "c", "a", "b", "c", "a", "b", "a", "b")
fitHigherOrder(sequence)
#> $lambda
#> [1] 1.000000e+00 1.852742e-09
#> 
#> $Q
#> $Q[[1]]
#>       a         b    c
#> a 0.125 0.4285714 0.75
#> b 0.750 0.1428571 0.25
#> c 0.125 0.4285714 0.00
#> 
#> $Q[[2]]
#>           a         b    c
#> a 0.1428571 0.5714286 0.25
#> b 0.4285714 0.2857143 0.75
#> c 0.4285714 0.1428571 0.00
#> 
#> 
#> $X
#>   a   b   c 
#> 0.4 0.4 0.2 
#> 
# weights by maximum likelihood (no Rsolnp needed)
fit <- fitHigherOrder(sequence, order = 2, method = "mle")
fit$lambda
#> [1] 1.000000e+00 2.350277e-08
higherOrderLogLik(sequence, fit)$logLik
#> [1] -13.08457
```
