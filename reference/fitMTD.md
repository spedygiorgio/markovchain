# Fit a mixture transition distribution (MTD) model

Estimates by maximum likelihood the mixture transition distribution
model of Raftery (1985) for a high-order Markov chain: the probability
of the next state is a weighted mixture of the contributions of the last
`order` states, all governed by the same transition matrix.

## Usage

``` r
fitMTD(
  sequence,
  order = 2,
  start = NULL,
  nstart = 1,
  tol = 1e-10,
  maxit = 10000L
)
```

## Arguments

- sequence:

  An empirical sequence of states (a character vector, or a vector
  coercible to character), without missing values.

- order:

  Order of the model, a positive integer smaller than the length of the
  sequence.

- start:

  Index of the first observation entering the likelihood, at least
  `order + 1` (the default).

- nstart:

  `nstart - 1` is the number of random starting points of the EM
  algorithm added to the deterministic ones (see Details).

- tol:

  Convergence tolerance on the relative change of the log-likelihood
  between two iterations.

- maxit:

  Maximum number of EM iterations for each starting point.

## Value

A list with components

- lambda:

  the estimated lag weights, named `lag1`, `lag2`, ...

- estimate:

  a `markovchain` object (by row) with the estimated transition matrix
  \\Q\\

- Q:

  a list of `order` copies of \\Q\\ stored by column
  (`Q[[g]][to, from]`), the layout returned by
  [`fitHigherOrder`](fitHigherOrder.md), so that the fit can be passed
  to [`higherOrderLogLik`](higherOrderLogLik.md)

- X:

  the relative frequencies of the states in `sequence`

- logLikelihood, npar, AIC, BIC, nobs:

  maximized log-likelihood, number of free parameters, information
  criteria and number of observations entering the likelihood

- order, start:

  as used in the fit

- iterations, converged:

  EM iterations and convergence flag of the returned fit

- model:

  the string `"MTD"`

## Details

For an order \\k\\ the model is \$\$P(X_t = j \mid X\_{t-1} = i_1,
\dots, X\_{t-k} = i_k) = \sum\_{g=1}^{k} \lambda_g\\ q\_{i_g j},\$\$
where \\Q = (q\_{ij})\\ is an \\r \times r\\ transition matrix (rows are
departure states) and the lag weights \\\lambda_g\\ sum to one. It needs
\\r(r-1) + k - 1\\ parameters instead of the \\r^k (r-1)\\ of a fully
parameterized Markov chain of order \\k\\, which makes high orders
practicable (Raftery, 1985; Berchtold and Raftery, 2002). This is
different from the model fitted by
[`fitHigherOrder`](fitHigherOrder.md), which mixes a different empirical
matrix for each lag and does not estimate them jointly with the weights.

The weights are constrained to be non-negative, as in most applications
of the model; Raftery's original formulation also allows negative
weights provided all the transition probabilities stay in \\\[0, 1\]\\,
which is not supported here. Under this constraint the likelihood is
maximized by the EM algorithm of Lebre and Bourguignon (2008), in which
the latent variable is the lag that generated each observation; it is
implemented in C++. Every iteration increases the likelihood, but the
likelihood of the MTD model can have several local maxima (Berchtold,
2001). The models of order \\1, \dots,\\ `order` are therefore fitted in
turn on the same observations, and order \\k\\ is started from equal
weights with the transition matrix of all lags pooled, and from the fit
of order \\k - 1\\ extended with a zero and with a small positive weight
for the new lag. The first of these extensions has the likelihood of
order \\k - 1\\, so the likelihood returned never decreases with the
order, as it must for nested models. `nstart - 1` further random
starting points can be added for the requested order; they are drawn at
random, so call [`set.seed`](https://rdrr.io/r/base/Random.html)
beforehand for reproducible results. The fit with the highest likelihood
is returned.

The likelihood is conditional on the observations before `start`: it is
the product of the transition probabilities of \\x_t\\ for \\t = \\
`start`, \\\dots\\, \\n\\. The default, `order + 1`, uses every
observation an order-`order` model can predict. To compare models of
different orders (or with the Markov chains of
[`fitHigherOrder`](fitHigherOrder.md) via
[`higherOrderLogLik`](higherOrderLogLik.md)) use the same `start` for
all, for instance \\1 +\\ the largest order; Berchtold and Raftery
(2002) condition on the first 14 observations, i.e. `start = 15`.

A state that never occurs in a conditioning position (for instance one
observed only at the end of the sequence) does not enter the likelihood;
its row of \\Q\\ is not identified and is returned as uniform.

AIC and BIC use the number of free parameters \\r(r-1) + k - 1\\, and
BIC the number of observations entering the likelihood. Berchtold and
Raftery (2002) do not count the elements of \\Q\\ estimated as exactly
zero; their BIC can be obtained by subtracting the number of those
elements from `npar`.

## References

Raftery, A. E. (1985). A model for high-order Markov chains. Journal of
the Royal Statistical Society, Series B, 47(3), 528-539.

Berchtold, A. (2001). Estimation in the mixture transition distribution
model. Journal of Time Series Analysis, 22(4), 379-397.

Berchtold, A. and Raftery, A. E. (2002). The mixture transition
distribution model for high-order Markov chains and non-Gaussian time
series. Statistical Science, 17(3), 328-356.

Lebre, S. and Bourguignon, P.-Y. (2008). An EM algorithm for estimation
in the mixture transition distribution model. Journal of Statistical
Computation and Simulation, 78(1), 1-15.

## See also

[`fitHigherOrder`](fitHigherOrder.md),
[`higherOrderLogLik`](higherOrderLogLik.md),
[`higherOrderPredict`](higherOrderPredict.md),
[`markovchainFit`](markovchainFit.md)

## Examples

``` r
# hourly wind directions at Koeberg (Berchtold and Raftery, 2002)
wind <- read.csv(system.file("extdata", "koeberg_wind.csv",
                             package = "markovchain"))$state
fit <- fitMTD(wind, order = 2, start = 15)
fit$lambda
#>      lag1      lag2 
#> 0.7568342 0.2431658 
fit$estimate
#> MTD(2) 
#>  A  4 - dimensional discrete Markov Chain defined by the following states: 
#>  1, 2, 3, 4 
#>  The transition matrix  (by rows)  is defined as follows: 
#>            1          2          3            4
#> 1 0.83067002 0.06867021 0.00767778 9.298199e-02
#> 2 0.03699661 0.90114666 0.06185673 7.502067e-55
#> 3 0.01535363 0.15522884 0.80712064 2.229689e-02
#> 4 0.07788932 0.00000000 0.05271010 8.694006e-01
#> 
c(logLik = fit$logLikelihood, BIC = fit$BIC)
#>    logLik       BIC 
#> -393.3968  872.5032 

# several starting points guard against local maxima
set.seed(1)
fitMTD(wind, order = 3, start = 15, nstart = 5)$logLikelihood
#> [1] -393.2398
```
