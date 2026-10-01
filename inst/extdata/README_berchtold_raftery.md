# Data from Berchtold and Raftery (2002)

Two categorical series used in Tables 2 and 3 of

Berchtold, A. and Raftery, A. E. (2002). The Mixture Transition Distribution
Model for High-Order Markov Chains and Non-Gaussian Time Series.
*Statistical Science*, 17(3), 328-356.

They are used by the package's unit tests and by the `higher_order_markov_chains`
and `an_introduction_to_markovchain_package` vignettes to check
`higherOrderLogLik()`, `assessIndependence()` and `fitMTD()` against the
published log-likelihoods and BIC values.

The files were kindly provided by the authors (Adrian E. Raftery and Andre
Berchtold) for use in this package; they stated that the data are not copyright
protected.

| file | rows | coding |
|------|------|--------|
| `koeberg_wind.csv` | 744 | `state`, hourly wind direction at Koeberg (South Africa) recoded into 4 directions, states `1`-`4`. The original series is from MacDonald and Zucchini (1997). |
| `epileptic_seizures.csv` | 204 | `seizure`, daily binary series: `0` = no epileptic seizure, `1` = at least one seizure. |

The authors' original epilepsy file is coded `1`/`2`; it was recoded to `0`/`1`
here (`1` -> `0`, `2` -> `1`), as they indicated.

## Convention of the published tables

Berchtold and Raftery condition the log-likelihoods on the first 14
observations, so that every model (independence, Markov chains of order up to
3, MTD models) is evaluated on the same `n - 14` observations; the BIC uses
`n - 14` as sample size, and the number of parameters counts only those that
are not forced to zero. Reproducing the tables therefore requires
`start = 15` in `higherOrderLogLik()` and `fitMTD()`.

MacDonald, I. L. and Zucchini, W. (1997). *Hidden Markov and Other Models for
Discrete-valued Time Series*. Chapman & Hall.
