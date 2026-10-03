# Select the order of a Markov chain by information criteria

Fits by maximum likelihood fully parameterized Markov chains of order
\\0, 1, \dots,\\ `maxOrder` to an empirical sequence (or to a list of
independent sequences), all on the same observations, and selects the
order that minimizes the BIC or the AIC (Tong, 1975; Katz, 1981).

## Usage

``` r
selectOrder(
  sequence,
  maxOrder = 3,
  criterion = c("BIC", "AIC"),
  start = NULL,
  parameters = c("full", "observed")
)
```

## Arguments

- sequence:

  An empirical sequence of states (a vector coercible to character,
  without missing values), or a list of such sequences.

- maxOrder:

  The highest order considered, a non-negative integer.

- criterion:

  The criterion used to select the order: `"BIC"` (default) or `"AIC"`.

- start:

  Index of the first observation of each sequence entering the
  likelihood, at least `maxOrder + 1` (the default).

- parameters:

  How the free parameters are counted: `"full"` (default) or
  `"observed"`, see Details.

## Value

A list with components

- order:

  the order selected by `criterion`

- criterion:

  the criterion used

- table:

  a data frame with one row per order: `order`, `logLik`, `npar`, `AIC`,
  `BIC`, and the likelihood-ratio test against the previous order (`LR`,
  `df`, `p.value`; `NA` for order 0)

- nobs:

  number of observations entering the likelihood

- start, parameters:

  as used

## Details

A Markov chain of order \\k\\ on \\r\\ states has a transition
probability for each context (the \\k\\ previous states) and next state;
order 0 is independence. Its maximum likelihood estimates are the
observed transition frequencies, and its log-likelihood is \\\sum N(c,
j) \log\\N(c, j) / N(c)\\\\, where \\N(c, j)\\ counts context \\c\\
followed by state \\j\\.

Information criteria are only comparable when every model is evaluated
on the same observations. An order-\\k\\ chain can only predict an
observation from position \\k + 1\\ onwards, so all orders are evaluated
on the observations from `start` (by default `maxOrder + 1`) to the end
of each sequence, the earlier ones being only used as contexts.
Berchtold and Raftery (2002), for instance, condition on the first 14
observations (`start = 15`). For a list of sequences the counts are
pooled, the observations before `start` of each sequence are only used
as contexts, and no context crosses from one sequence to the next.

With `parameters = "full"` (the default) an order-\\k\\ chain has \\r^k
(r - 1)\\ free parameters, as in Tong (1975) and Katz (1981). With
`parameters = "observed"` only the probabilities that are not estimated
as zero are counted, that is the number of distinct next states minus
one summed over the observed contexts; this is the convention of
Berchtold and Raftery (2002), which is less penalizing when many
transitions are never observed. In both cases BIC uses the number of
observations entering the likelihood.

The table also reports the likelihood-ratio statistic of each order
against the previous one, \\G = 2 (\ell_k - \ell\_{k-1})\\, with degrees
of freedom equal to the difference in the number of parameters and an
asymptotic chi-squared p-value (Anderson and Goodman, 1957). These tests
are not adjusted for multiplicity, and they and the criteria become
unreliable when the number of contexts \\r^k\\ is not small compared
with the number of observations: BIC is consistent for the order
(Csiszar and Shields, 2000), whereas AIC tends to select too high an
order in long sequences (Katz, 1981).

## References

Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about
Markov chains. The Annals of Mathematical Statistics, 28(1), 89-110.

Tong, H. (1975). Determination of the order of a Markov chain by
Akaike's information criterion. Journal of Applied Probability, 12(3),
488-497.

Katz, R. W. (1981). On some criteria for estimating the order of a
Markov chain. Technometrics, 23(3), 243-249.

Csiszar, I. and Shields, P. C. (2000). The consistency of the BIC Markov
order estimator. The Annals of Statistics, 28(6), 1601-1619.

Berchtold, A. and Raftery, A. E. (2002). The mixture transition
distribution model for high-order Markov chains and non-Gaussian time
series. Statistical Science, 17(3), 328-356.

## See also

[`assessOrder`](statisticalTests.md),
[`verifyMarkovProperty`](statisticalTests.md),
[`fitHigherOrder`](fitHigherOrder.md), [`fitMTD`](fitMTD.md),
[`higherOrderLogLik`](higherOrderLogLik.md)

## Examples

``` r
# Alofi rainfall: three states
data(rain)
selectOrder(rain$rain, maxOrder = 3)$table
#>   order    logLik npar      AIC      BIC        LR df      p.value
#> 1     0 -1133.827    2 2271.655 2281.648        NA NA           NA
#> 2     1 -1038.063    6 2088.125 2118.105 191.52944  4 2.486724e-40
#> 3     2 -1025.095   18 2086.190 2176.130  25.93541 12 1.096202e-02
#> 4     3 -1005.563   54 2119.127 2388.948  39.06299 36 3.338344e-01

# Koeberg wind directions with the conventions of Berchtold and Raftery (2002)
wind <- read.csv(system.file("extdata", "koeberg_wind.csv",
                             package = "markovchain"))$state
sel <- selectOrder(wind, maxOrder = 3, start = 15, parameters = "observed")
sel$order
#> [1] 1
sel$table
#>   order    logLik npar       AIC       BIC         LR df       p.value
#> 1     0 -954.8297    3 1915.6594 1929.4385         NA NA            NA
#> 2     1 -413.3036   11  848.6073  899.1308 1083.05209  8 1.751228e-228
#> 3     2 -374.9209   27  803.8418  927.8540   76.76548 16  6.334196e-10
#> 4     3 -346.1945   39  770.3890  949.5178   57.45277 12  6.547231e-08
```
