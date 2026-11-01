# Statistical tests for Markov chains

Tests of hypotheses on observed sequences of states, based on transition
counts (Anderson and Goodman 1957; Kullback, Kupperman and Ku 1962):

- `verifyMarkovProperty()`: the first-order Markov property, i.e. that
  the next state is independent of the previous one given the present
  one (against a second-order chain);

- `assessOrder()`: the same hypothesis with the Pearson statistic, also
  for a list of independent sequences;

- `verifyEmpiricalToTheoretical()`: that the sequence (or a matrix of
  transition counts) comes from a given chain;

- `verifyHomogeneity()`: that several sequences (or count matrices)
  share the same transition matrix;

- `assessStationarity()`: that the transition probabilities do not
  change over time, comparing `nblocks` consecutive blocks of the
  sequence.

## Usage

``` r
verifyMarkovProperty(
  sequence,
  method = c("G", "Pearson", "simulation"),
  B = 9999,
  seed = NULL,
  verbose = TRUE
)

assessOrder(sequence, verbose = TRUE)

verifyEmpiricalToTheoretical(
  data,
  object,
  method = c("G", "Pearson", "simulation"),
  B = 9999,
  seed = NULL,
  verbose = TRUE
)

verifyHomogeneity(
  inputList,
  method = c("G", "Pearson", "simulation"),
  B = 9999,
  seed = NULL,
  verbose = TRUE
)

assessStationarity(sequence, nblocks, structural.zeros = NULL, verbose = TRUE)
```

## Arguments

- sequence:

  An observed sequence of states (a vector); for `assessOrder()` also a
  list of independent sequences, whose triples of consecutive states are
  formed within each sequence and pooled.

- method:

  The statistic: `"G"` (likelihood ratio, the default), `"Pearson"`, or
  `"simulation"` (likelihood ratio with a Monte Carlo p-value).

- B:

  Number of Monte Carlo samples for `method = "simulation"`.

- seed:

  Optional random seed for `method = "simulation"`.

- verbose:

  Logical: print the result?

- data:

  An observed sequence of states, or a matrix of transition counts with
  the states as dimnames.

- object:

  A [`markovchain`](markovchain-class.md) object, the chain of the null
  hypothesis.

- inputList:

  A list of at least two observed sequences or matrices of transition
  counts (with the states as dimnames).

- nblocks:

  Number of consecutive blocks of (nearly) equal length into which the
  sequence is split, at least two.

- structural.zeros:

  Optional logical matrix, with one row and column per state (in the
  order of `unique(sequence)`, or with the states as dimnames), marking
  the transitions that are impossible; they are left out of the
  comparison.

## Value

An object of class `"htest"`, printed if `verbose = TRUE` and returned
invisibly, with the usual components `statistic`, `parameter` (degrees
of freedom), `p.value`, `method` and `data.name`, and also `dof`.
Depending on the test it also contains the `observed` and `expected`
counts and the estimated `transitionMatrix` (`verifyMarkovProperty()`,
`verifyEmpiricalToTheoretical()`), the `pooled` counts and the
`transitionCounts` of each element (`verifyHomogeneity()`), and, with
`method = "simulation"`, the simulated statistics (`simulations`) and
`B`.

## Details

Each test compares observed transition counts with the counts expected
under the null hypothesis, by the likelihood-ratio statistic \\G^2 = 2
\sum O \log(O / E)\\ (`method = "G"`) or by Pearson's \\X^2 = \sum (O -
E)^2 / E\\ (`method = "Pearson"`); cells with no expected count are left
out, and the degrees of freedom count only the transitions actually
observed. The p-value is from the chi-squared distribution, or, with
`method = "simulation"`, from `B` samples simulated under the null
hypothesis (the \\G^2\\ statistic, p-value \\(1 + \\\\simulated \ge
observed\\) / (B + 1)\\), which is more reliable when many counts are
small. Sequences of numbers or factors are compared as character
strings.

## References

Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about
Markov chains. The Annals of Mathematical Statistics 28(1), 89-110.

Kullback, S., Kupperman, M. and Ku, H. H. (1962). Tests for contingency
tables and Markov chains. Technometrics 4(4), 573-608.

## See also

[`assessIndependence()`](assessIndependence.md) (independence of
consecutive states), [`selectOrder()`](selectOrder.md) (order by
information criteria), [`compareChains()`](compareChains.md).

## Examples

``` r
data(rain, package = "markovchain")
# first-order Markov property of the Alofi rainfall sequence
verifyMarkovProperty(rain$rain)
#> 
#>  Likelihood-ratio test
#> 
#> data:  rain$rain
#> G-squared = 25.837, df = 12, p-value = 0.01132
#> 
assessOrder(rain$rain)
#> 
#>  Pearson's Chi-squared test for Markov order
#> 
#> data:  rain$rain
#> X-squared = 26.096, df = 12, p-value = 0.0104
#> 

# compatibility with a chain: the fitted one, and a different one
fit <- markovchainFit(rain$rain)$estimate
verifyEmpiricalToTheoretical(rain$rain, fit)
#> 
#>  Likelihood-ratio test
#> 
#> data:  rain$rain
#> G-squared = 0, df = 6, p-value = 1
#> 
s <- c("0", "1-5", "6+")
uniform <- new("markovchain", states = s,
               transitionMatrix = matrix(1 / 3, 3, 3, dimnames = list(s, s)))
verifyEmpiricalToTheoretical(rain$rain, uniform)
#> 
#>  Likelihood-ratio test
#> 
#> data:  rain$rain
#> G-squared = 325.12, df = 6, p-value < 2.2e-16
#> 

# the first and second half of the sequence have the same transitions
verifyHomogeneity(list(rain$rain[1:548], rain$rain[549:1096]))
#> 
#>  Likelihood-ratio test
#> 
#> data:  list(rain$rain[1:548], rain$rain[549:1096])
#> G-squared = 11.231, df = 6, p-value = 0.0815
#> 
# time-homogeneity over four blocks
assessStationarity(rain$rain, nblocks = 4)
#> 
#>  Pearson's Chi-squared test for time-homogeneity
#> 
#> data:  rain$rain
#> X-squared = 29.051, df = 18, p-value = 0.04776
#> 

# a Monte Carlo p-value
verifyMarkovProperty(rain$rain, method = "simulation", B = 499, seed = 1)
#> 
#>  Parametric Monte Carlo test (likelihood-ratio statistic)
#> 
#> data:  rain$rain
#> G-squared = 25.837, df = 12, p-value = 0.008
#> 
```
