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
  maxit = 10000L,
  model = c("MTD", "MTDg"),
  accelerate = FALSE,
  se = FALSE,
  negativeWeights = FALSE
)
```

## Arguments

- sequence:

  An empirical sequence of states (a character vector, or a vector
  coercible to character), without missing values, or a list of such
  vectors, taken as independent sequences (see Details).

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

  Maximum number of EM iterations for each starting point (EM steps, if
  `accelerate = TRUE`).

- model:

  `"MTD"` (default), the model with a single matrix, or `"MTDg"`, the
  general model with a different matrix for each lag.

- accelerate:

  Logical: speed up the EM algorithm with the SQUAREM scheme of Varadhan
  and Roland (2008)? The default `FALSE` keeps the plain algorithm.
  Convergence is much faster when the likelihood is flat, and the same
  maximum is reached.

- se:

  Logical: compute standard errors (see Details)?

- negativeWeights:

  Logical: allow negative lag weights, as in Raftery (1985)? Only for
  `model = "MTD"`; the default `FALSE` keeps them non-negative.

## Value

A list with components

- lambda:

  the estimated lag weights, named `lag1`, `lag2`, ...

- estimate:

  a `markovchain` object (by row) with the estimated transition matrix
  \\Q\\ (for `model = "MTDg"`, a list with the `order` objects of the
  matrices \\Q_1, \dots, Q_k\\)

- Q:

  a list of `order` copies of \\Q\\ (the \\Q_g\\ for `"MTDg"`) stored by
  column (`Q[[g]][to, from]`), the layout returned by
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

  EM iterations (sweeps of the coordinate ascent with
  `negativeWeights = TRUE`) and convergence flag of the returned fit

- model:

  `"MTD"` or `"MTDg"`

- negativeWeights:

  as used in the fit

- se:

  only if `se = TRUE`: a list with the standard errors `lambda` of the
  weights, `Q` of the matrix (same layout as `estimate`; `NA` for the
  entries that cannot be estimated and, for MTDg, for all the weights
  and matrices) and `probabilities` of the transition probabilities, a
  matrix with one row per history observed in the data (named by its
  states, oldest first, separated by `"-"`, as passed to
  [`higherOrderPredict`](higherOrderPredict.md)) and one column per next
  state; see Details

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

By default the weights are constrained to be non-negative, as in most
applications of the model (see below for negative weights). Under this
constraint the likelihood is maximized by the EM algorithm of Lebre and
Bourguignon (2008), in which the latent variable is the lag that
generated each observation; it is implemented in C++. Every iteration
increases the likelihood, but the likelihood of the MTD model can have
several local maxima (Berchtold, 2001). The models of order \\1,
\dots,\\ `order` are therefore fitted in turn on the same observations,
and order \\k\\ is started from equal weights with the transition matrix
of all lags pooled, and from the fit of order \\k - 1\\ extended with a
zero and with a small positive weight for the new lag. The first of
these extensions has the likelihood of order \\k - 1\\, so the
likelihood returned never decreases with the order, as it must for
nested models. `nstart - 1` further random starting points can be added
for the requested order; they are drawn at random, so call
[`set.seed`](https://rdrr.io/r/base/Random.html) beforehand for
reproducible results. The fit with the highest likelihood is returned.

The likelihood is conditional on the observations before `start`: it is
the product of the transition probabilities of \\x_t\\ for \\t = \\
`start`, \\\dots\\, \\n\\. The default, `order + 1`, uses every
observation an order-`order` model can predict. To compare models of
different orders (or with the Markov chains of
[`fitHigherOrder`](fitHigherOrder.md) via
[`higherOrderLogLik`](higherOrderLogLik.md)) use the same `start` for
all, for instance \\1 +\\ the largest order; Berchtold and Raftery
(2002) condition on the first 14 observations, i.e. `start = 15`.

With a list of sequences the likelihood is the product over the
sequences: each contributes the observations from `start` to its end,
the earlier ones being only used as conditioning states, no lag spans
two sequences, and a sequence with fewer than `start` states is not used
(apart from the state frequencies `X`). `order` must be smaller than the
length of the longest sequence.

A state that never occurs in a conditioning position (for instance one
observed only at the end of the sequence) does not enter the likelihood;
its row of \\Q\\ is not identified and is returned as uniform.

A warning is given when the number of observations entering the
likelihood does not exceed the number of free parameters.

With `model = "MTDg"` (Berchtold, 2001; Lebre and Bourguignon, 2008) the
probability of the next state is \\\sum_g \lambda_g\\ q^{(g)}\_{i_g j}\\
with a transition matrix \\Q_g\\ for each lag. It has \\(r-1)(1 +
k(r-1))\\ free parameters, the number of the models of
[`fitHigherOrder`](fitHigherOrder.md), of which it is the maximum
likelihood version: there the matrices are the empirical ones and only
the weights are fitted. The EM algorithm is the same; order \\k\\ is
started from the MTD fit of the same order, from the empirical matrices
with the best weights for them, and from the fit of order \\k - 1\\, so
its likelihood is never below those of `model = "MTD"` and of the
`method = "mle"` fit of [`fitHigherOrder`](fitHigherOrder.md), nor of
lower orders. The row of \\Q_g\\ of a state never seen at lag \\g\\ is
returned as uniform.

The weights and matrices of the MTDg model are not identified: since the
probabilities only depend on the products \\A_g = \lambda_g Q_g\\,
adding \\\epsilon\\ times the same probability vector to every row of
\\A_1\\ and subtracting it from every row of \\A_2\\ changes
\\\lambda_1\\, \\\lambda_2\\, \\Q_1\\ and \\Q_2\\ but no transition
probability. The maximum is therefore attained on a set of dimension up
to \\(k - 1) r\\, of which the values returned are one point; only the
transition probabilities (and the likelihood) are determined, and `npar`
counts the identified parameters. The lag weights of an MTDg fit should
not be interpreted as the importance of the lags.

With `negativeWeights = TRUE` (model `"MTD"` only) the weights may be
negative, as in Raftery's (1985) original formulation, provided every
transition probability \\\sum_g \lambda_g q\_{i_g j}\\ lies in \\\[0,
1\]\\ for all the combinations of departure states \\i_g\\: for each
arrival state \\j\\, \\\min_i q\_{ij}\\ times the sum of the positive
weights must not be below \\\max_i q\_{ij}\\ times the sum of the
absolute negative weights (Raftery, 1985). The EM algorithm does not
apply, and the likelihood is maximized by block coordinate ascent in the
spirit of Berchtold (2001): in the weights and in each row of \\Q\\ in
turn, mass is moved from the entry with the smallest to the one with the
largest partial derivative, by the step that maximizes the likelihood
along that direction within the feasible region. Order \\k\\ starts from
the fit with non-negative weights of the same order and from the fit of
order \\k - 1\\ with a zero weight for the new lag, so its likelihood is
never below either of them; for the requested order it also starts with
each lag in turn given a negative weight. The likelihood has several
local maxima, and this is slower than the EM algorithm. `iterations`
then counts the sweeps of the ascent at the requested order.
Unidentified rows of \\Q\\ are returned as uniform, or as the mean of
the other rows if a uniform row would violate the condition.

With `se = TRUE` the standard errors are those of the observed
information at the maximum, by the delta method, obtained in closed form
with the entries of each probability vector (the weights and each row of
\\Q\\) summing to one. Entries that are on the boundary of the parameter
space (below \\10^{-4}\\) are held fixed, as in the usual treatment of
constrained maxima, and have `NA`; a warning is given if a weight is
among them, since then the standard errors of the others are conditional
on that weight being zero (with negative weights allowed, if Raftery's
condition is binding). Standard errors are also given for the transition
probabilities of every history observed in the data. For the MTDg model
only these are available, the weights and matrices not being identified;
they are computed with a generalized inverse of the information, which
gives the same variance for any identified quantity. All are asymptotic,
and approximate in small samples or with sparse \\Q\\.

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

Varadhan, R. and Roland, C. (2008). Simple and globally convergent
methods for accelerating the convergence of any EM algorithm.
Scandinavian Journal of Statistics, 35(2), 335-353.

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

# standard errors, and the model with a matrix for each lag
fit2 <- fitMTD(wind, order = 2, start = 15, se = TRUE)
fit2$se$lambda
#>       lag1       lag2 
#> 0.05006497 0.05006497 
fitMTD(wind, order = 2, start = 15, model = "MTDg")$logLikelihood
#> [1] -389.4452
```
