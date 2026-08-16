# Published literature benchmarks for Markov-chain tests

## A published non-rejection example for the Markov property: Burkett, Scherer and Todd (2019)

Burkett, Scherer and Todd study recurrent and persistent states in
financial markets and use chi-squared tests to investigate the
first-order Markov property ([Burkett et al.
2019](#ref-burkett2019recurrent)). Their published examples are useful
benchmarks because the observed counts used in the tests are reported
explicitly.

For one of their tests, the transition of interest is
$\left. 1\rightarrow 1 \right.$. Observations are grouped according to
the state immediately preceding the origin state. The published counts
are:

``` r
burkett <- matrix(
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

burkett
```

                     transition to state 1 transition elsewhere
    previous state 1                   127                   59
    previous state 2                    11                    7
    previous state 4                    36                   20
    previous state 5                     9                    4

The Pearson test on these published counts does not reject the null
hypothesis that the probability of the transition
$\left. 1\rightarrow 1 \right.$ is independent of the preceding state.
The Pearson calculation gives $\chi^{2} \approx 0.64$ with 3 degrees of
freedom and a p-value of about 0.89, comfortably above 0.05. This
provides an independent, published example in which a Markov-property
hypothesis is **not rejected**.

``` r
burkett.test <- stats::chisq.test(burkett, correct = FALSE)
```

    Warning in stats::chisq.test(burkett, correct = FALSE): Chi-squared
    approximation may be incorrect

``` r
burkett.test
```

        Pearson's Chi-squared test

    data:  burkett
    X-squared = 0.63611, df = 3, p-value = 0.8881

This example is deliberately presented as a literature benchmark for the
contingency-table test underlying Markov-property inference. It should
not be interpreted as a direct call to
[`verifyMarkovProperty()`](../reference/statisticalTests.md) on a
sequence: the paper performs a cell-wise conditional test, whereas
[`verifyMarkovProperty()`](../reference/statisticalTests.md) evaluates
the first-order Markov property jointly over the observed
triple-transition table.

## Published empirical-versus-theoretical transition matrix: Kullback, Kupperman and Ku (1962)

Kullback, Kupperman and Ku ([Kullback et al.
1962](#ref-kullback1962tests)) develop likelihood-ratio and
information-theoretic tests for finite Markov chains, including the
hypothesis that an empirical transition matrix is compatible with a
**specified theoretical transition matrix**. Their paper contains a
worked example that is reproduced in the historical `markovchain`
documentation and provides exactly the type of benchmark needed for
[`verifyEmpiricalToTheoretical()`](../reference/statisticalTests.md).

The published sequence is:

``` r
kullback.sequence <- c(
  0,1,2,2,1,0,0,0,0,0,0,1,2,2,2,1,0,0,1,0,0,0,0,0,0,1,1,
  2,0,0,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,2,1,0,
  0,2,1,0,0,0,0,0,0,1,1,1,2,2,0,0,2,1,1,1,1,2,1,1,1,1,1,1,1,1,1,0,2,
  0,1,1,0,0,0,0,0,0,0,0,2,2,2,0,0,0,1,1,1,1,1,2,1,2,0,0,0,1,1,
  0,0,0,0,0,2,2,1,1,1,1,1,2,1,2,0,0,0,1,2,2,2,0,0,0,1,1
)

kullback.counts <- createSequenceMatrix(kullback.sequence)
kullback.counts
```

       0  1  2
    0 54 11  7
    1 10 29 10
    2  7 10  9

The theoretical transition matrix specified in the published example is

``` r
kullback.P <- matrix(
  c(5/8, 1/4, 1/8,
    1/4, 1/2, 1/4,
    1/4, 3/8, 3/8),
  byrow = TRUE,
  nrow = 3,
  dimnames = list(as.character(0:2), as.character(0:2))
)

theoretical <- new(
  "markovchain",
  states = as.character(0:2),
  transitionMatrix = kullback.P
)

theoretical
```

    Unnamed Markov chain
     A  3 - dimensional discrete Markov Chain defined by the following states:
     0, 1, 2
     The transition matrix  (by rows)  is defined as follows:
          0     1     2
    0 0.625 0.250 0.125
    1 0.250 0.500 0.250
    2 0.250 0.375 0.375

The published calculation gives a likelihood-ratio chi-squared statistic
of approximately $6.55$ on 6 degrees of freedom, with p-value
approximately **0.3643**. Hence the null hypothesis that the empirical
transition matrix is compatible with the specified theoretical matrix is
**not rejected** at the 5% level ([Kullback et al.
1962](#ref-kullback1962tests)).

The same calculation can be reproduced with the package:

``` r
result <- verifyEmpiricalToTheoretical(
  data = kullback.sequence,
  object = theoretical,
  method = "G",
  verbose = TRUE
)
```

        Likelihood-ratio test

    data:  kullback.sequence
    G-squared = 7.1034, df = 6, p-value = 0.3114

``` r
result
```

        Likelihood-ratio test

    data:  kullback.sequence
    G-squared = 7.1034, df = 6, p-value = 0.3114

The historical `markovchain` article reproduces the transition counts as
$\begin{pmatrix}
51 & 11 & 8 \\
12 & 31 & 9 \\
6 & 11 & 10
\end{pmatrix}$ and reports the likelihood-ratio statistic 6.551795, 6
degrees of freedom and p-value 0.3642899 ([Kullback et al.
1962](#ref-kullback1962tests)). This gives us an unusually useful
regression benchmark: the empirical transition matrix is not rejected as
compatible with the specified theoretical Markov chain.

For comparison, using the Pearson version gives the same inferential
conclusion:

``` r
verifyEmpiricalToTheoretical(
  data = kullback.sequence,
  object = theoretical,
  method = "Pearson",
  verbose = TRUE
)
```

        Pearson's Chi-squared test

    data:  kullback.sequence
    X-squared = 6.7223, df = 6, p-value = 0.3473

This example is preferable to a simulated example when demonstrating
[`verifyEmpiricalToTheoretical()`](../reference/statisticalTests.md):
both the empirical sequence and the theoretical transition matrix are
taken from a published worked example, and the null hypothesis is
explicitly **not rejected**.

## Why these examples are useful together

The examples provide complementary published benchmarks for the
package’s statistical-inference functionality:

- **Burkett et al. (2019):** an independently published, real-data
  example in which a Markov-property hypothesis is not rejected.
- **Kullback et al. (1962):** a published empirical-versus-theoretical
  transition-matrix test with $p \approx 0.364$, so the specified Markov
  chain is not rejected.
- **Tran and Carmichael (2012):** a real contractor-payment application
  in the main statistical-analysis vignette in which time-homogeneity is
  rejected.

Together they avoid presenting only rejection examples and demonstrate
the important inferential principle that a large p-value means **failure
to reject the null hypothesis**, not proof that the Markov model is
true.

Burkett, Matthew W., William T. Scherer, and Andrew Todd. 2019. “Are
Financial Market States Recurrent and Persistent?” *Cogent Economics and
Finance* 7 (1): 1622171.
<https://doi.org/10.1080/23322039.2019.1622171>.

Kullback, Solomon, Michael Kupperman, and H. H. Ku. 1962. “Tests for
Contingency Tables and Markov Chains.” *Technometrics* 4 (4): 573–608.
<https://doi.org/10.1080/00401706.1962.10490041>.
