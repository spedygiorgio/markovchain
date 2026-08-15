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
verifyMarkovProperty(rain$rain, method = "G", verbose = FALSE)
verifyMarkovProperty(rain$rain, method = "Pearson", verbose = FALSE)
```

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
verifyEmpiricalToTheoretical(x, mc, method = "G", verbose = FALSE)
```

The function also accepts a raw transition-count matrix:

``` r
counts <- createSequenceMatrix(x)
verifyEmpiricalToTheoretical(counts, mc, method = "Pearson", verbose = FALSE)
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
  seed = 123, verbose = FALSE
)$p.value
#> [1] 0.29
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
verifyHomogeneity(list(x1, x2, x3), method = "G", verbose = FALSE)
```

Pearson and likelihood-ratio versions are available:

``` r
verifyHomogeneity(list(x1, x2, x3), method = "Pearson", verbose = FALSE)
```

A parametric Monte Carlo version simulates multinomial rows under the
pooled transition matrix while retaining the observed number of
transitions originating from each state in each sample. The pooled model
is re-estimated for each replicate:

``` r
verifyHomogeneity(
  list(x1, x2, x3), method = "simulation", B = 999,
  seed = 123, verbose = FALSE
)$p.value
#> [1] 0.213
```

These tests are useful when transition data come from different periods,
regions, portfolios, or cohorts and one wants to assess whether a common
homogeneous transition matrix is adequate. The likelihood-ratio
formulation follows the classical contingency-table and Markov-chain
framework of Kullback et al. ([1962](#ref-kullback1962tests)).

### Published example: Alofi daily rainfall

Avery and Henderson ([P. J. Avery and D. A. Henderson
1999](#ref-averyHenderson)) analysed discrete-state time series
including the Alofi daily rainfall data distributed with `markovchain`.
The data classify successive days into three rainfall states and are
included as `rain`.

``` r
data(rain)
table(rain$rain)
#> 
#>   0 1-5  6+ 
#> 548 295 253
createSequenceMatrix(rain$rain)
#>       0 1-5  6+
#> 0   362 126  60
#> 1-5 136  90  68
#> 6+   50  79 124
```

A published likelihood comparison of first- and second-order Markov
models for these data reports maximized log-likelihoods of approximately
-1038.06 and -1025.10, respectively ([Davison
2003](#ref-davison2003statistical)). Thus

$$G^{2} = 2\{(-1025.10) - (-1038.06)\} = 25.92.$$

For three states the comparison has 12 degrees of freedom, giving an
asymptotic p-value of approximately 0.011. Thus the first-order model is
rejected at the 5% level in that likelihood comparison.

The package performs the corresponding first-order goodness-of-fit
calculation directly from the sequence:

``` r
verifyMarkovProperty(rain$rain, method = "G", verbose = FALSE)
```

The exact value need not equal the rounded published likelihood
calculation, but the hypothesis is the same: whether the state
immediately preceding the current rainfall state contains additional
information about the next state.

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
                     seed = 123, verbose = FALSE)$p.value
#> [1] 0.018
verifyEmpiricalToTheoretical(counts, mc, method = "simulation", B = 999,
                             seed = 123, verbose = FALSE)$p.value
#> [1] 0.29
verifyHomogeneity(list(x1, x2, x3), method = "simulation", B = 999,
                  seed = 123, verbose = FALSE)$p.value
#> [1] 0.213
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

Davison, A. C. 2003. *Statistical Models*. Cambridge University Press.
<https://doi.org/10.1017/CBO9780511815850>.

Kullback, S, Michael Kupperman, and HH Ku. 1962. “Tests for Contingency
Tables and Marltov Chains.” *Technometrics* 4 (4): 573–608.

P. J. Avery, and D. A. Henderson. 1999. “Fitting Markov Chain Models to
Discrete State Series.” *Applied Statistics* 48 (1): 53–61.
