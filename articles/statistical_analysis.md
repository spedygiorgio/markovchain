# Statistical analysis of Markov chain sequences

## Statistical analysis

### Testing the first-order Markov property

For an observed sequence $X_{1},\ldots,X_{n}$ with $r$ possible states,
the first-order Markov property requires that

$$P(X_{t + 1} = k \mid X_{t} = j,X_{t - 1} = i) = P(X_{t + 1} = k \mid X_{t} = j).$$

Equivalently, $X_{t - 1}$ and $X_{t + 1}$ are conditionally independent
given $X_{t}$. Let $N_{ijk}$ denote the number of consecutive triples
$(i,j,k)$. Under the first-order null model,

$$E_{ijk} = N_{ij +}\frac{N_{+jk}}{N_{+j +}}.$$

This formulation is closely related to the classical likelihood-based
inference for Markov chains ([Anderson and Goodman
1957](#ref-anderson1957statistical); [Kullback et al.
1962](#ref-kullback1962tests)). The package implements the Pearson
statistic

$$X^{2} = \sum\limits_{i,j,k}\frac{(N_{ijk} - E_{ijk})^{2}}{E_{ijk}}$$

and the likelihood-ratio statistic

$$G^{2} = 2\sum\limits_{N_{ijk} > 0}N_{ijk}{\log}\left( \frac{N_{ijk}}{E_{ijk}} \right).$$

For complete support the asymptotic degrees of freedom are
$r(r - 1)^{2}$. With sparse support the implementation derives the
degrees of freedom from the observed conditional support. A
non-rejection means only that the data provide insufficient evidence
against the first-order model; it does not prove the Markov property.

``` r
data(rain)
verifyMarkovProperty(rain$rain, method = "G", verbose = TRUE)
#> 
#>  Likelihood-ratio test
#> 
#> data:  rain$rain
#> G-squared = 25.837, df = 12, p-value = 0.01132
verifyMarkovProperty(rain$rain, method = "Pearson", verbose = TRUE)
#> 
#>  Pearson's Chi-squared test
#> 
#> data:  rain$rain
#> X-squared = 26.096, df = 12, p-value = 0.0104
```

The printed output above reports the statistic, degrees of freedom and
p-value for the two asymptotic tests. Keeping `verbose = TRUE` in the
vignette makes the numerical result part of the rendered example rather
than silently suppressing it.

### Empirical versus theoretical transition probabilities

[`verifyEmpiricalToTheoretical()`](../reference/statisticalTests.md)
tests whether observed transition counts are compatible with a specified
transition matrix. For each origin state $i$, the null model is
multinomial with probabilities given by row $i$ of the theoretical
transition matrix. The likelihood-ratio and Pearson versions use the
same observed row totals and differ only in their goodness-of-fit
statistic. This is the natural row-wise multinomial formulation of
goodness-of-fit for transition probabilities ([Kullback et al.
1962](#ref-kullback1962tests)).

``` r
P <- matrix(c(.70, .20, .10,
              .20, .50, .30,
              .10, .30, .60), 3, 3, byrow = TRUE,
            dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
mc <- new("markovchain", states = c("a", "b", "c"), transitionMatrix = P)
set.seed(1)
x <- rmarkovchain(1000, mc, t0 = "a")
verifyEmpiricalToTheoretical(x, mc, method = "G", verbose = TRUE)
#> 
#>  Likelihood-ratio test
#> 
#> data:  x
#> G-squared = 7.2259, df = 6, p-value = 0.3005
```

The function also accepts a raw transition-count matrix:

``` r
counts <- createSequenceMatrix(x)
verifyEmpiricalToTheoretical(counts, mc, method = "Pearson", verbose = TRUE)
#> 
#>  Pearson's Chi-squared test
#> 
#> data:  counts
#> X-squared = 7.3534, df = 6, p-value = 0.2894
```

If an observed transition has positive count while the theoretical
probability is zero, the theoretical model is incompatible with the
data. Conversely, a zero observed count where the theoretical
probability is positive is not by itself evidence against the model.

For small or sparse samples, the asymptotic chi-squared approximation
can be replaced by a parametric Monte Carlo calibration. The observed
row totals are retained and multinomial transition counts are simulated
under the theoretical matrix:

``` r
verifyEmpiricalToTheoretical(
  counts, mc, method = "simulation", B = 999,
  seed = 123, verbose = TRUE
)
#> 
#>  Parametric Monte Carlo test (likelihood-ratio statistic)
#> 
#> data:  counts
#> G-squared = 7.2259, df = 6, p-value = 0.29
```

The resulting p-value uses the usual Monte Carlo correction

$$\frac{1 + B_{\geq}}{B + 1}.$$

This is a parametric Monte Carlo test, not an exact finite-sample test.

### Homogeneity across samples

[`verifyHomogeneity()`](../reference/statisticalTests.md) tests whether
several empirical sequences or transition-count matrices share the same
transition probabilities. Under the null hypothesis, a common transition
matrix is estimated by pooling the transition counts across samples. The
test can therefore be viewed, row by row, as a multinomial homogeneity
test.

``` r
set.seed(123)
x1 <- rmarkovchain(500, mc, t0 = "a")
x2 <- rmarkovchain(500, mc, t0 = "b")
x3 <- rmarkovchain(500, mc, t0 = "c")
verifyHomogeneity(list(x1, x2, x3), method = "G", verbose = TRUE)
#> 
#>  Likelihood-ratio test
#> 
#> data:  list(x1, x2, x3)
#> G-squared = 15.475, df = 12, p-value = 0.2165
```

Pearson and likelihood-ratio versions are available:

``` r
verifyHomogeneity(list(x1, x2, x3), method = "Pearson", verbose = TRUE)
#> 
#>  Pearson's Chi-squared test
#> 
#> data:  list(x1, x2, x3)
#> X-squared = 14.21, df = 12, p-value = 0.2875
```

A parametric Monte Carlo version simulates multinomial rows under the
pooled transition matrix while retaining the observed number of
transitions originating from each state in each sample. The pooled model
is re-estimated for each replicate:

``` r
verifyHomogeneity(
  list(x1, x2, x3), method = "simulation", B = 999,
  seed = 123, verbose = TRUE
)
#> 
#>  Parametric Monte Carlo test (likelihood-ratio statistic)
#> 
#> data:  list(x1, x2, x3)
#> G-squared = 15.475, df = 12, p-value = 0.213
```

These tests are useful when transition data come from different periods,
regions, portfolios, or cohorts and one wants to assess whether a common
homogeneous transition matrix is adequate. The likelihood-ratio
formulation follows the classical contingency-table and Markov-chain
framework of Kullback et al. ([1962](#ref-kullback1962tests)).

### Testing time-homogeneity of transition probabilities

[`assessStationarity()`](../reference/statisticalTests.md) tests a
different form of stability: whether the transition probabilities are
constant across a specified number of consecutive blocks of the
sequence. For each origin state, the transition counts in the blocks
form a contingency table. The null hypothesis is that the conditional
distribution of the next state, given the current state, is the same in
every informative block.

This distinction is important.
[`verifyHomogeneity()`](../reference/statisticalTests.md) compares
**different samples** and asks whether they share a common transition
matrix, whereas
[`assessStationarity()`](../reference/statisticalTests.md) splits **one
sequence into consecutive blocks** and asks whether its transition
matrix is stable over time. The latter is therefore a test of
time-homogeneity of a first-order Markov chain.

The test is based on the Pearson chi-squared statistic. Transition
counts are retained as counts; zero counts caused merely by finite
sampling are not treated as impossible transitions. Blocks in which a
particular origin state is never observed contain no information about
that state’s conditional distribution and are therefore omitted from
that state’s contingency table.

``` r
set.seed(123)
stationary.sequence <- rmarkovchain(1000, mc, t0 = "a")
assessStationarity(stationary.sequence, nblocks = 4, verbose = TRUE)
#> 
#>  Pearson's Chi-squared test for time-homogeneity
#> 
#> data:  stationary.sequence
#> X-squared = 15.05, df = 18, p-value = 0.6585
```

The output gives the combined Pearson statistic, its degrees of freedom,
and the asymptotic p-value. A large p-value does not establish
stationarity; it means that the data do not provide sufficient evidence
against constant transition probabilities across the blocks.

#### Structural zeros

A zero transition count can have two very different meanings. A
**sampling zero** is simply a transition that did not occur in the
observed sample. A **structural zero** is a transition that is
impossible by construction of the process. Only the latter should be
excluded from the contingency table and from the corresponding degrees
of freedom.

[`assessStationarity()`](../reference/statisticalTests.md) therefore
accepts the optional `structural.zeros` argument. It is a logical
transition matrix with `TRUE` marking transitions that are structurally
impossible. Its rows and columns correspond to the states of the
sequence; when dimnames are supplied they are used to match the states
explicitly.

For example, suppose state `a` can never be followed by `a`:

``` r
structural.zeros <- matrix(FALSE, 3, 3,
                           dimnames = list(c("a", "b", "c"),
                                           c("a", "b", "c")))
structural.zeros["a", "a"] <- TRUE

structural.sequence <- c("a", "b", "a", "c", "a", "b", "c", "a",
                         "b", "c", "b", "a", "c", "b", "a", "c",
                         "a", "b", "c", "a", "b", "c", "b", "a")

assessStationarity(
  structural.sequence,
  nblocks = 4,
  structural.zeros = structural.zeros,
  verbose = TRUE
)
#> 
#>  Pearson's Chi-squared test for time-homogeneity
#> 
#> data:  structural.sequence
#> X-squared = 1.7639, df = 15, p-value = 1
```

A structural-zero transition that is actually observed is rejected as
inconsistent with the supplied model:

``` r
invalid.zeros <- matrix(FALSE, 3, 3,
                        dimnames = list(c("a", "b", "c"),
                                        c("a", "b", "c")))
invalid.zeros["a", "b"] <- TRUE

try(
  assessStationarity(
    structural.sequence,
    nblocks = 4,
    structural.zeros = invalid.zeros,
    verbose = TRUE
  )
)
#> Error in .assessStationarity(sequence, nblocks = nblocks, structural.zeros = structural.zeros,  : 
#>   sequence contains a transition marked as a structural zero
```

The distinction between sampling and structural zeros is especially
relevant for sparse transition matrices. Automatically interpreting
every observed zero as a structural zero would change the null model and
can substantially inflate or otherwise distort the degrees of freedom.

### Published real-data examples

The stationarity of transition probabilities has been tested on real
data for decades. Anderson and Goodman (1957) provide the classical
likelihood-ratio framework, while Collins (1974) applies Markov-chain
inference to manufacturing establishments and tests the homogeneity of
transition matrices observed over successive years ([Anderson and
Goodman 1957](#ref-anderson1957statistical); [Collins
1974](#ref-collins1974micro)). Tran and Carmichael (2012) apply the
likelihood-ratio approach specifically to contractor payment histories
and test the assumption that transition probabilities are constant over
project time ([Tran and Carmichael 2012](#ref-tranCarmichael2012)).

#### Manufacturing establishments: Collins (1974)

Collins ([Collins 1974](#ref-collins1974micro)) applies tests of the
Markov property, chain order and **homogeneity of transition matrices**
to data on manufacturing establishments. The study uses micro-unit data
from manufacturing establishments observed over successive years. The
published analysis finds no statistically significant evidence against
homogeneity for the annual transition matrices. This provides a useful
real-data example in which the null hypothesis of constant transition
probabilities is *not* rejected, complementing examples in which
stationarity is rejected.

The original article is the appropriate source for the historical data
and complete transition tables; the package vignette does not reproduce
proprietary or incompletely documented micro-unit records. We therefore
use the Collins result as a published benchmark rather than
manufacturing a pseudo-sequence from rounded tables.

#### Contractor payments: Tran and Carmichael (2012)

Tran and Carmichael ([Tran and Carmichael
2012](#ref-tranCarmichael2012)) provide two real construction-project
case studies. Their null hypothesis is

$$H_{0}:p_{jk}(\tau) = p_{jk},\qquad\tau = 1,\ldots,s,$$

against time-dependent transition probabilities. They use the
Anderson–Goodman likelihood-ratio test and note that the Pearson
chi-squared homogeneity test is asymptotically equivalent. Both case
studies reject stationarity at the 5% level. The road project has 41
progress claims over 32 months; the steel-fabrication project has 12
claims over 12 months ([Tran and Carmichael
2012](#ref-tranCarmichael2012)).

The second case is particularly useful for `markovchain` because the
paper publishes the empirical, month-by-month outstanding amounts. Table
4 reports the cumulative total claimed amount and amounts outstanding by
more than 2, 4, 6 and 8 weeks. The values below reproduce the published
table in thousands of dollars ([Tran and Carmichael
2012](#ref-tranCarmichael2012)).

``` r
steel <- data.frame(
  month = 1:12,
  total = c(89531, 173329, 421008, 531654, 1082009, 1358504,
            1657696, 1709917, 1801566, 1818082, 1847810, 1888942),
  over2 = c(89531, 173329, 198096, 308742, 344784, 407462,
            440712, 492933, 492933, 509449, 509449, 550581),
  over4 = c(89531, 99710, 124477, 134641, 170683, 233361,
            266611, 267075, 267075, 267075, 267075, 308207),
  over6 = c(89531, 99710, 124477, 134641, 170683, 233361,
            266611, 267075, 267075, 267075, 267075, 270101),
  over8 = c(89531, 99710, 124477, 134641, 170683, 233361,
            266611, 267075, 267075, 267075, 267075, 270101)
)

steel
#>    month   total  over2  over4  over6  over8
#> 1      1   89531  89531  89531  89531  89531
#> 2      2  173329 173329  99710  99710  99710
#> 3      3  421008 198096 124477 124477 124477
#> 4      4  531654 308742 134641 134641 134641
#> 5      5 1082009 344784 170683 170683 170683
#> 6      6 1358504 407462 233361 233361 233361
#> 7      7 1657696 440712 266611 266611 266611
#> 8      8 1709917 492933 267075 267075 267075
#> 9      9 1801566 492933 267075 267075 267075
#> 10    10 1818082 509449 267075 267075 267075
#> 11    11 1847810 509449 267075 267075 267075
#> 12    12 1888942 550581 308207 270101 270101
```

The paper converts these outstanding amounts into transition flows
$\alpha_{jk}(\tau)$ and then applies the likelihood-ratio statistic. Its
reported result is

$$-2{\log}(\lambda) = 1636.1,\qquad df = 6 \times 5 \times 11 = 330,$$

with a p-value below 0.05. The null hypothesis of time-homogeneous
transition probabilities is therefore rejected ([Tran and Carmichael
2012](#ref-tranCarmichael2012)).

Because these are **monetary transition flows rather than individual
claim-level state observations**,
[`assessStationarity()`](../reference/statisticalTests.md) should not be
applied to the twelve rows as if each dollar were an independent
observation. Instead, the same contingency-table calculation is
represented by [`verifyHomogeneity()`](../reference/statisticalTests.md)
using the published transition-flow matrices. The following code
reconstructs the six-state transition-count matrices from Table 4. The
two absorbing states are included as self-transitions only to complete
the transition matrices; they do not carry information about time
variation.

``` r
steel.transition <- function(z) {
  m <- matrix(0, 6, 6,
              dimnames = list(as.character(0:5), as.character(0:5)))

  # States 0--3 are progressively older outstanding states;
  # states 4 and 5 are Paid and To be resolved.
  m[1, 2] <- z$over2
  m[1, 5] <- z$total - z$over2
  m[2, 3] <- z$over4
  m[2, 5] <- z$over2 - z$over4
  m[3, 4] <- z$over6
  m[3, 5] <- z$over4 - z$over6
  m[4, 5] <- z$over6 - z$over8
  m[4, 6] <- z$over8

  # Absorbing states.
  m[5, 5] <- 1
  m[6, 6] <- 1
  m
}

steel.matrices <- lapply(seq_len(nrow(steel)), function(i) {
  steel.transition(steel[i, ])
})

verifyHomogeneity(
  steel.matrices,
  method = "G",
  verbose = TRUE
)
#> 
#>  Likelihood-ratio test
#> 
#> data:  steel.matrices
#> G-squared = 1186103, df = 33, p-value < 2.2e-16
```

The reconstructed-data calculation gives the same **qualitative
conclusion** as Tran and Carmichael: the p-value is overwhelmingly small
and the hypothesis of constant transition probabilities is rejected. The
numerical statistic is not expected to equal 1636.1 exactly, because the
paper’s likelihood-ratio construction uses the monetary-flow likelihood
and its own phase/time aggregation, whereas
[`verifyHomogeneity()`](../reference/statisticalTests.md) applies the
package’s general contingency-table formulation to the reconstructed
transition-flow matrices. The comparison is therefore deliberately made
at the level of the inferential conclusion rather than by asserting
numerical identity.

This example also illustrates why the new `structural.zeros` facility is
useful. In the owner-payment model, only transitions to the next ageing
state or to an absorbing state are admissible; many cells of the
$6 \times 6$ transition matrix are structural zeros ([Tran and
Carmichael 2012](#ref-tranCarmichael2012)). For a genuine
individual-level sequence from the same model, these constraints can be
supplied directly to
[`assessStationarity()`](../reference/statisticalTests.md) through a
logical transition matrix.

``` r
payment.structural.zeros <- matrix(TRUE, 6, 6,
  dimnames = list(as.character(0:5), as.character(0:5)))

# Allowed transitions: ageing by one state, payment, and resolution;
# states 4 and 5 are absorbing.
payment.structural.zeros[1, c(2, 5)] <- FALSE
payment.structural.zeros[2, c(3, 5)] <- FALSE
payment.structural.zeros[3, c(4, 5)] <- FALSE
payment.structural.zeros[4, c(5, 6)] <- FALSE
payment.structural.zeros[5, 5] <- FALSE
payment.structural.zeros[6, 6] <- FALSE

payment.structural.zeros
#>      0     1     2     3     4     5
#> 0 TRUE FALSE  TRUE  TRUE FALSE  TRUE
#> 1 TRUE  TRUE FALSE  TRUE FALSE  TRUE
#> 2 TRUE  TRUE  TRUE FALSE FALSE  TRUE
#> 3 TRUE  TRUE  TRUE  TRUE FALSE FALSE
#> 4 TRUE  TRUE  TRUE  TRUE FALSE  TRUE
#> 5 TRUE  TRUE  TRUE  TRUE  TRUE FALSE
```

The two published case studies therefore provide complementary real-data
evidence: Collins reports a case in which homogeneity is not rejected,
whereas Tran and Carmichael report two construction-payment cases in
which it is rejected. They also illustrate an important practical point:
structural zeros are part of the model specification and should not be
confused with transitions that merely happen not to appear in a finite
sample.

### Monte Carlo calibration and validation

The three statistical tests provide asymptotic G and Pearson tests and,
where appropriate, a parametric Monte Carlo alternative. Monte Carlo
calibration is particularly useful when expected counts are small or the
transition support is sparse.

For [`verifyMarkovProperty()`](../reference/statisticalTests.md), a
first-order transition matrix is fitted under the null, sequences are
simulated from it, and the first-order model is re-fitted in every
replicate. For
[`verifyEmpiricalToTheoretical()`](../reference/statisticalTests.md) and
[`verifyHomogeneity()`](../reference/statisticalTests.md), the
simulation conditions on the observed row totals and samples multinomial
transition counts under the theoretical or pooled transition
probabilities. In all cases the Monte Carlo p-value is

$$\frac{1 + \sum\limits_{b = 1}^{B}I(T_{b} \geq T_{obs})}{B + 1},$$

where $T$ is the selected test statistic.

For example:

``` r
verifyMarkovProperty(rain$rain, method = "simulation", B = 999,
                     seed = 123, verbose = TRUE)
#> 
#>  Parametric Monte Carlo test (likelihood-ratio statistic)
#> 
#> data:  rain$rain
#> G-squared = 25.837, df = 12, p-value = 0.018
verifyEmpiricalToTheoretical(counts, mc, method = "simulation", B = 999,
                             seed = 123, verbose = TRUE)
#> 
#>  Parametric Monte Carlo test (likelihood-ratio statistic)
#> 
#> data:  counts
#> G-squared = 7.2259, df = 6, p-value = 0.29
verifyHomogeneity(list(x1, x2, x3), method = "simulation", B = 999,
                  seed = 123, verbose = TRUE)
#> 
#>  Parametric Monte Carlo test (likelihood-ratio statistic)
#> 
#> data:  list(x1, x2, x3)
#> G-squared = 15.475, df = 12, p-value = 0.213
```

These procedures are parametric Monte Carlo tests, not exact
finite-sample tests. A separate validation study is deliberately kept
outside `R CMD check` because it is computationally expensive. It should
evaluate Type-I error, power, sample-size sensitivity, state-space size,
sparse transition matrices, and calibration of the asymptotic tests
against the Monte Carlo alternatives. The validation script is stored
under `inst/validation/` and can be run manually.

### References

Anderson, T. W., and L. A. Goodman. 1957. “Statistical Inference about
Markov Chains.” *The Annals of Mathematical Statistics* 28 (1): 89–110.
<https://doi.org/10.1214/aoms/1177707039>.

Collins, Lyndhurst. 1974. “Estimating Markov Transition Probabilities
from Micro-Unit Data.” *Journal of the Royal Statistical Society. Series
C (Applied Statistics)* 23 (3): 355–71.
<https://doi.org/10.2307/2347128>.

Kullback, S, Michael Kupperman, and HH Ku. 1962. “Tests for Contingency
Tables and Marltov Chains.” *Technometrics* 4 (4): 573–608.

Tran, Hanh, and David G. Carmichael. 2012. “Stationarity of the
Transition Probabilities in the Markov Chain Formulation of Owner
Payments on Projects.” *ANZIAM Journal* 53: C69–89.
<https://doi.org/10.21914/anziamj.v53i0.5111>.
