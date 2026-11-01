# Fit a higher order Markov chain

Fits the higher order Markov chain of Ching et al., in which the
distribution of the next state is a weighted mixture of the empirical
lag-\\i\\ transition matrices, \\i = 1, \dots,\\ `order`.
`seq2freqProb()` returns the state frequencies of a sequence and
`seq2matHigh()` its empirical lag-`order` transition matrix, the
building blocks of the model.

## Usage

``` r
fitHigherOrder(sequence, order = 2, method = c("mle", "lsq"))
seq2freqProb(sequence)
seq2matHigh(sequence, order)
```

## Arguments

- sequence:

  The sequence of states: a character vector, or a vector coercible to
  character (numbers and factors are matched to the states as character
  strings), without missing values. A list of such vectors is taken as
  independent sequences, see Details.

- order:

  Order of the chain, a positive integer smaller than the sequence
  length; for `seq2matHigh()`, the lag.

- method:

  How the weights \\\lambda\\ are estimated: `"mle"` (default) or
  `"lsq"`, see Details.

## Value

`fitHigherOrder()`: a list with `lambda`, the weights of the lags; `Q`,
the list of the lag-\\i\\ transition matrices (column-stochastic: entry
\\\[j, k\]\\ is the probability of moving to \\j\\ from \\k\\); and `X`,
the state frequencies. With `method = "lsq"`, `NULL` with a warning if
Rsolnp is not available.

`seq2freqProb()`: the vector of state frequencies, named by state.
`seq2matHigh()`: the column-stochastic lag-`order` transition matrix (a
state never observed as a departure state gets a uniform column).

## Details

The fitted model expresses the distribution of the next state as the
mixture \\\sum\_{i=1}^{k} \lambda_i Q_i x\_{t-i}\\ of the empirical
lag-\\i\\ transition matrices \\Q_i\\ (see `seq2matHigh`), with weights
\\\lambda_i \ge 0\\ summing to one. The matrices \\Q_i\\ are the same
for both methods; only the weights differ.

`method = "lsq"` (the default before version 1.2.1, and the only method
before the argument existed) chooses \\\lambda\\ to minimize the squared
distance between the stationary distribution and its image under the
mixture, as in Ching et al.; it needs the Rsolnp package and returns
`NULL` with a warning if it is unavailable. This criterion is weak: each
lag matrix maps the empirical distribution onto itself up to end
effects, \\Q_i X \approx X\\ with an error of order \\i/n\\, so the
objective is nearly flat in \\\lambda\\ and the weights it returns can
be unstable from one sample to another. It remains available, to
reproduce earlier results or the method of Ching et al., but it is no
longer the default.

`method = "mle"` (the default) chooses \\\lambda\\ to maximize the
log-likelihood \\\sum\_{t=k+1}^{n} \log \sum_i \lambda_i Q_i\[x_t,
x\_{t-i}\]\\ of the observations a model of order \\k\\ can predict. For
fixed \\Q_i\\ the problem is concave, so the maximum is global, and it
is solved, without Rsolnp, by Newton's method with an active set on the
simplex, which stops when the Karush-Kuhn-Tucker conditions of the
maximum hold within \\10^{-10}\\ (in practice after 3 to 10 iterations,
also when some weights are 0). The weights are therefore those that give
the highest value of [`higherOrderLogLik`](higherOrderLogLik.md) for the
same observations. (Up to version 1.2 the EM algorithm for mixture
weights was used: near a weight equal to 0, or when the lags carry
nearly the same information, it converged very slowly and could stop
visibly short of the maximum.) They are maximum likelihood estimates
*conditional on* the empirical matrices \\Q_i\\, which are not
re-estimated: this is not the maximum likelihood estimator of the
mixture with free matrices, in which the weights are in general not
identifiable (a common distribution can be moved from \\\lambda_j Q_j\\
to \\\lambda_i Q_i\\ without changing any transition probability), so
the weights should not be read as the relative importance of the lags
beyond this conditional sense. Note that this is not the mixture
transition distribution model of Raftery (1985), in which a single
matrix is shared by all lags and is estimated together with the weights;
that model is fitted by [`fitMTD`](fitMTD.md).

With a list of sequences the lag matrices \\Q_i\\ and the state
frequencies are computed from the counts pooled over the sequences,
which are treated as independent: no pair of observations spans two
sequences, and a sequence with no more than `order` states contributes
to the frequencies only. The weights are then estimated as for a single
sequence, from the observations that a model of that order can predict
in every sequence.

A warning is given when the sequence supplies no more observations than
the model has free parameters (the number used by
[`higherOrderLogLik`](higherOrderLogLik.md)), since the weights are then
essentially undetermined.

## References

Ching, W. K., Huang, X., Ng, M. K., & Siu, T. K. (2013). Higher-order
markov chains. In Markov Chains (pp. 141-176). Springer US.

Ching, W.-K., Ng, M. K. and Fung, E. S. (2008). Higher-order
multivariate Markov chains and their applications. Linear Algebra and
its Applications 428(2-3), 492-507.
[doi:10.1016/j.laa.2007.05.021](https://doi.org/10.1016/j.laa.2007.05.021)

Raftery, A. E. (1985). A model for high-order Markov chains. Journal of
the Royal Statistical Society, Series B, 47(3), 528-539.

## See also

[`higherOrderLogLik()`](higherOrderLogLik.md),
[`higherOrderPredict()`](higherOrderPredict.md),
[`fitMTD()`](fitMTD.md), [`selectOrder()`](selectOrder.md)

Other higher order chains:
[`fitHighOrderMultivarMC()`](fitHighOrderMultivarMC.md),
[`predictHommc()`](predictHommc.md)

## Author

Giorgio Spedicato, Tae Seung Kang

## Examples

``` r
sequence<-c("a", "a", "b", "b", "a", "c", "b", "a", "b", "c", "a", "b",
            "c", "a", "b", "c", "a", "b", "a", "b")
# weights by maximum likelihood (the default, no Rsolnp needed)
fit <- fitHigherOrder(sequence, order = 2)
fit$lambda
#> [1] 1 0
# weights by least squares on the stationary distribution, as in Ching et al.
fitHigherOrder(sequence, order = 2, method = "lsq")
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
higherOrderLogLik(sequence, fit)$logLik
#> [1] -13.08457
```
