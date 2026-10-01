# Test independence of consecutive states of an empirical sequence

Tests the null hypothesis that consecutive observations are independent,
\\P(X\_{t+1} = j \mid X_t = i) = P(X\_{t+1} = j)\\ for all \\i, j\\,
against the alternative of a first-order Markov chain. This is the
classical test of Anderson and Goodman (1957): a chi-squared (or
likelihood-ratio) test of independence applied to the table of observed
one-step transition counts, whose rows are the state at time \\t\\ and
whose columns are the state at time \\t + 1\\.

## Usage

``` r
assessIndependence(sequence, method = c("Pearson", "G"), verbose = TRUE)
```

## Arguments

- sequence:

  An empirical sequence of states (at least three observations, without
  missing values).

- method:

  Test statistic: `"Pearson"` (default) for the chi-squared statistic or
  `"G"` for the likelihood-ratio statistic.

- verbose:

  Should test results be printed?

## Value

An `htest` object, returned invisibly, with the additional components
`observed` (transition counts, rows are departure states) and `expected`
(counts expected under independence).

## Details

Only states actually observed as a departure state (rows) or as an
arrival state (columns) contribute to the degrees of freedom, which are
\\(r - 1)(c - 1)\\ for \\r\\ such rows and \\c\\ such columns. When
\\r\\ or \\c\\ is 1 (for instance a constant sequence) the test is not
defined: the degrees of freedom are 0 and the p-value is `NA`.

The statistic is asymptotic and, as for any chi-squared test on a table
of counts, unreliable when many expected counts are small (say below 5).
Successive transitions overlap (each observation is the arrival state of
one transition and the departure state of the next); this is the
standard treatment of the Anderson-Goodman test and is asymptotically
valid under the null hypothesis.

## References

Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about
Markov chains. *The Annals of Mathematical Statistics*, 28(1), 89–110.

## See also

[`verifyMarkovProperty`](statisticalTests.md),
[`assessOrder`](statisticalTests.md),
[`assessStationarity`](statisticalTests.md)

## Examples

``` r
# an independent sequence: the test should not reject
set.seed(1)
iid <- sample(c("a", "b", "c"), 500, replace = TRUE)
assessIndependence(iid)
#> 
#>  Pearson's Chi-squared test for independence of consecutive states
#> 
#> data:  iid
#> X-squared = 8.9139, df = 4, p-value = 0.06329
#> 

# a strongly dependent (Markov) sequence: the test rejects
mc <- new("markovchain", states = c("a", "b"),
          transitionMatrix = matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE))
dep <- rmarkovchain(500, mc, t0 = "a")
assessIndependence(dep, method = "G")
#> 
#>  Likelihood-ratio test for independence of consecutive states
#> 
#> data:  dep
#> G-squared = 245.08, df = 1, p-value < 2.2e-16
#> 
```
