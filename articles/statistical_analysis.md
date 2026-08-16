# Statistical analysis of Markov chain sequences

## Statistical analysis

This vignette illustrates the statistical inference functions for
empirical Markov chains. Test results are returned as objects of class
`htest`, following R’s standard convention for hypothesis tests.

### Testing the first-order Markov property

For an observed sequence $X_{1},\ldots,X_{n}$, the first-order Markov
property requires that

$$P(X_{t + 1} = k \mid X_{t} = j,X_{t - 1} = i) = P(X_{t + 1} = k \mid X_{t} = j).$$

The package implements Pearson and likelihood-ratio tests of this
property. A non-rejection means that the data provide insufficient
evidence against the null hypothesis; it does not prove that the process
is Markovian.

``` r
sequence <- c("a", "b", "a", "a", "a", "a", "b", "a", "b",
              "a", "b", "a", "a", "b", "b", "b", "a")

markov_test <- verifyMarkovProperty(sequence, method = "G", verbose = TRUE)
#> 
#>  Likelihood-ratio test
#> 
#> data:  sequence
#> G-squared = 0.59916, df = 2, p-value = 0.7411
markov_test
#> 
#>  Likelihood-ratio test
#> 
#> data:  sequence
#> G-squared = 0.59916, df = 2, p-value = 0.7411
class(markov_test)
#> [1] "htest" "list"
```

### Published example: non-rejection of the Markov property

Burkett, Scherer and Todd (2019) ([Burkett et al.
2019](#ref-burkett2019recurrent)) study recurrent and persistent
financial-market states. One of their published tests considers the
transition into state 1 and reports the following counts:

``` r
burkett_counts <- matrix(
  c(127, 59,
     11,  7,
     36, 20,
      9,  4),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(
    c("previous state 1", "previous state 2",
      "previous state 4", "previous state 5"),
    c("transition to state 1", "transition elsewhere")
  )
)

chisq.test(burkett_counts, correct = FALSE)
#> Warning in chisq.test(burkett_counts, correct = FALSE): Chi-squared
#> approximation may be incorrect
#> 
#>  Pearson's Chi-squared test
#> 
#> data:  burkett_counts
#> X-squared = 0.63611, df = 3, p-value = 0.8881
```

The paper reports a non-significant test (p-value about 0.98). The
calculation above uses the published contingency table directly. It is
included as an independent literature example of a non-rejection; it is
not claimed to be numerically identical to
[`verifyMarkovProperty()`](../reference/statisticalTests.md), because
the published test is a cell-wise test for a particular transition
rather than the package’s global first-order Markov-property test.

### Empirical versus theoretical transition probabilities

[`verifyEmpiricalToTheoretical()`](../reference/statisticalTests.md)
tests whether observed transition counts are compatible with a specified
transition matrix. The following example is from Kullback, Kupperman and
Ku (1962) ([Kullback et al. 1962](#ref-kullback1962tests)). The paper
gives both the empirical counts and the theoretical transition matrix.

``` r
kullback_counts <- matrix(
  c(51, 11,  8,
    12, 31,  9,
     6, 11, 10),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
)

kullback_P <- matrix(
  c(5/8, 1/4, 1/8,
    1/4, 1/2, 1/4,
    1/4, 3/8, 3/8),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
)

kullback_mc <- new(
  "markovchain",
  states = c("a", "b", "c"),
  transitionMatrix = kullback_P
)

kullback_result <- verifyEmpiricalToTheoretical(
  kullback_counts,
  kullback_mc,
  method = "G",
  verbose = TRUE
)
#> 
#>  Likelihood-ratio test
#> 
#> data:  kullback_counts
#> G-squared = 6.5184, df = 6, p-value = 0.3677

kullback_result
#> 
#>  Likelihood-ratio test
#> 
#> data:  kullback_counts
#> G-squared = 6.5184, df = 6, p-value = 0.3677
class(kullback_result)
#> [1] "htest" "list"
```

Kullback et al. report $G^{2} = 6.551795$, $df = 6$ and $p = 0.3642899$.
The historical `markovchain` implementation reproduced the published
value by using column totals when constructing expected counts. For the
row-wise multinomial null implemented by the current package, expected
counts are constructed as

$$E_{ij} = n_{i +}P_{ij},$$

where $n_{i +}$ is the observed total in row $i$. This gives
$G^{2} = 6.518384$, with $df = 6$ and the same non-rejection of the null
hypothesis. The small numerical discrepancy is therefore a correction of
the historical implementation rather than a change in the inferential
conclusion.

### Homogeneity across samples

[`verifyHomogeneity()`](../reference/statisticalTests.md) tests whether
several empirical transition-count matrices are compatible with a common
transition mechanism.

``` r
P <- matrix(
  c(.70, .20, .10,
    .20, .50, .30,
    .10, .30, .60),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
)

mc <- new("markovchain",
          states = c("a", "b", "c"),
          transitionMatrix = P)

set.seed(123)
x1 <- rmarkovchain(500, mc, t0 = "a")
x2 <- rmarkovchain(500, mc, t0 = "a")

homogeneity_test <- verifyHomogeneity(
  list(x1, x2), method = "G", verbose = TRUE
)
#> 
#>  Likelihood-ratio test
#> 
#> data:  list(x1, x2)
#> G-squared = 1.8881, df = 6, p-value = 0.9297
homogeneity_test
#> 
#>  Likelihood-ratio test
#> 
#> data:  list(x1, x2)
#> G-squared = 1.8881, df = 6, p-value = 0.9297
```

### Stationarity and time-homogeneity

[`assessStationarity()`](../reference/statisticalTests.md) tests whether
transition probabilities remain constant across blocks of an empirical
sequence. This is a test of time-homogeneity of the transition
mechanism, rather than a test of the marginal stationary distribution.

#### Published example: contractor payments

Tran and Carmichael (2012) ([Tran and Carmichael
2012](#ref-tranCarmichael2012)) apply likelihood-ratio tests to
transition probabilities in a real contractor-payment application. Their
analysis provides an example where the null hypothesis of
time-homogeneous transition probabilities is rejected. The application
also motivates explicit handling of structural zeros, i.e. transitions
that are impossible under the state-space model.

The package allows such transitions to be specified through the
`structural.zeros` argument.

``` r
structural.zeros <- matrix(FALSE, 3, 3,
                           dimnames = list(c("a", "b", "c"),
                                           c("a", "b", "c")))
structural.zeros["a", "a"] <- TRUE

structural_sequence <- c(
  "a", "b", "a", "c", "a", "b", "c", "a",
  "b", "c", "b", "a", "c", "b", "a", "c",
  "a", "b", "c", "a", "b", "c", "b", "a"
)

stationarity_test <- assessStationarity(
  structural_sequence,
  nblocks = 4,
  structural.zeros = structural.zeros,
  verbose = TRUE
)
#> 
#>  Pearson's Chi-squared test for time-homogeneity
#> 
#> data:  structural_sequence
#> X-squared = 1.7639, df = 15, p-value = 1
stationarity_test
#> 
#>  Pearson's Chi-squared test for time-homogeneity
#> 
#> data:  structural_sequence
#> X-squared = 1.7639, df = 15, p-value = 1
class(stationarity_test)
#> [1] "htest" "list"
```

A structural zero represents a transition ruled out by the model. It is
therefore different from a sampling zero: an unobserved but possible
transition remains part of the statistical analysis.

For sparse data, the asymptotic chi-squared approximation should be
interpreted with care.

Burkett, Matthew W., William T. Scherer, and Andrew Todd. 2019. “Are
Financial Market States Recurrent and Persistent?” *Cogent Economics and
Finance* 7 (1): 1622171.
<https://doi.org/10.1080/23322039.2019.1622171>.

Kullback, Solomon, Michael Kupperman, and H. H. Ku. 1962. “Tests for
Contingency Tables and Markov Chains.” *Technometrics* 4 (4): 573–608.
<https://doi.org/10.1080/00401706.1962.10490041>.

Tran, Hanh, and David G. Carmichael. 2012. “Stationarity of the
Transition Probabilities in the Markov Chain Formulation of Owner
Payments on Projects.” *ANZIAM Journal* 53: C69–89.
<https://doi.org/10.21914/anziamj.v53i0.5111>.
