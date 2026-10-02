# Log-likelihood, deviance and information criteria of a higher order Markov chain

Evaluates the log-likelihood of an empirical sequence under the higher
order Markov chain returned by [`fitHigherOrder`](fitHigherOrder.md),
and derives the deviance, AIC and BIC, so that models of different
orders can be compared.

## Usage

``` r
higherOrderLogLik(sequence, fit = NULL, order = 2, start = NULL)
```

## Arguments

- sequence:

  The empirical sequence of states, a character vector or a vector
  coercible to character (numbers and factors are matched to the states
  of the fit as character strings).

- fit:

  The list returned by [`fitHigherOrder`](fitHigherOrder.md) (or by
  [`fitMTD`](fitMTD.md)) for `sequence`. If `NULL`,
  `fitHigherOrder(sequence, order)` is computed.

- order:

  Order of the model to fit when `fit` is `NULL` (ignored otherwise; the
  order is then `length(fit$lambda)`).

- start:

  Index of the first observation included in the likelihood. Defaults to
  `order + 1`, the earliest observation a model of that order can
  predict.

## Value

A list with components `logLik`, `deviance`, `AIC`, `BIC`, `nobs`
(number of observations entering the likelihood), `npar`, `order` and
`start`. The log-likelihood is `-Inf` if an observed transition has
probability zero under the model (possible only for a sequence other
than the one the model was fitted on).

## Details

The fitted model is the mixture-transition-distribution model of Raftery
(1985) in the form used by Ching et al.: the probability of moving to
state \\x_t\\ given the past is \$\$P(x_t \mid x\_{t-1}, \dots,
x\_{t-k}) = \sum\_{i=1}^{k} \lambda_i\\ Q_i\[x_t, x\_{t-i}\],\$\$ where
\\Q_i\\ is the lag-\\i\\ transition matrix (see
[`seq2matHigh`](fitHigherOrder.md)) and \\k\\ is the order. The
log-likelihood is the sum of the logarithms of these probabilities over
the observations \\t = \\ `start`, \\\dots\\, \\T\\, and the deviance is
\\-2\\ times the log-likelihood.

Two points matter when interpreting the output. First, `fitHigherOrder`
chooses \\\lambda\\ by default (`method = "lsq"`) by least squares on
the stationary distribution, not by maximum likelihood, so the value
returned is the log-likelihood *of the fitted model*, not the maximum
attainable one; with `method = "mle"` the weights maximize this
log-likelihood for the observations a model of that order can predict.
Second, a model of order \\k\\ can only be evaluated from observation
\\k + 1\\ onwards; to compare orders on exactly the same data set
`start` to \\1 +\\ the largest order compared, otherwise the models are
evaluated on different numbers of observations and neither the
log-likelihood nor the information criteria are comparable.

The number of parameters used for AIC and BIC is the dimension of the
set of transition laws the model can represent, \\(r - 1)(1 + k (r -
1))\\, with \\r\\ the number of states: each of the \\k\\ lag matrices
has \\r (r - 1)\\ free probabilities and there are \\k - 1\\ free
weights, but the weights of a mixture of lag matrices are not
identifiable (a distribution common to all departure states can be moved
from \\\lambda_j Q_j\\ to \\\lambda_i Q_i\\ without changing any
transition probability), which removes \\r (k - 1)\\ parameters from the
naive count \\k\\ r (r - 1) + (k - 1)\\: for each next state the
transition probability is a sum of one term for each lag, a main-effects
function of the \\k\\ past states, with \\1 + k (r - 1)\\ free
coefficients. For \\k = 1\\ both counts are \\r (r - 1)\\. For a fit
returned by [`fitMTD`](fitMTD.md), whose lags share a single matrix, it
is \\r (r - 1) + (k - 1)\\.

## References

Raftery, A. E. (1985). A model for high-order Markov chains. Journal of
the Royal Statistical Society, Series B, 47(3), 528-539.

Ching, W. K., Huang, X., Ng, M. K., & Siu, T. K. (2013). Higher-order
markov chains. In Markov Chains (pp. 141-176). Springer US.

## See also

[`fitHigherOrder`](fitHigherOrder.md), [`fitMTD`](fitMTD.md),
[`higherOrderPredict`](higherOrderPredict.md)

## Examples

``` r
sequence <- c("a", "a", "b", "b", "a", "c", "b", "a", "b", "c", "a", "b",
              "c", "a", "b", "c", "a", "b", "a", "b")
# compare orders 1 and 2 on the same observations (start = 3)
if (requireNamespace("Rsolnp", quietly = TRUE)) {
  fit1 <- fitHigherOrder(sequence, order = 1)
  fit2 <- fitHigherOrder(sequence, order = 2)
  sapply(list(order1 = fit1, order2 = fit2), function(f)
    unlist(higherOrderLogLik(sequence, f, start = 3)[c("logLik", "deviance", "AIC", "BIC")]))
}
#>             order1    order2
#> logLik   -13.08457 -13.08457
#> deviance  26.16914  26.16914
#> AIC       38.16914  46.16914
#> BIC       43.51137  55.07286
```
