# Test the first-order Markov property of an empirical sequence

Tests the null hypothesis that the conditional distribution of the next
state depends only on the current state, against a second-order
alternative.

Tests whether the sequence is compatible with a first-order Markov chain
against a second-order alternative, by testing independence of past and
future states conditional on the present state. Degrees of freedom are
summed, present-state by present-state, over only the past and future
states actually observed with that present state, mirroring the other
functions documented on this page.

Tests whether transition probabilities are constant across consecutive
blocks. Structural zeros can be supplied explicitly through a logical
transition matrix.

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

  An empirical sequence.

- method:

  Test statistic: \`"G"\`, \`"Pearson"\`, or \`"simulation"\`.

- B:

  Number of Monte Carlo replicates for simulation.

- seed:

  Optional random seed.

- verbose:

  Should test results be printed?

- data:

  An empirical sequence or a matrix of transition counts.

- object:

  A \`markovchain\` object specifying theoretical probabilities.

- inputList:

  A list whose elements are empirical sequences or matrices.

- nblocks:

  Number of blocks, at least two.

- structural.zeros:

  Optional logical matrix marking impossible transitions.

## Value

An \`htest\` object with additional package-specific components.

An \`htest\` object.

An \`htest\` object with observed and expected counts.

An \`htest\` object with pooled and individual transition counts.

An \`htest\` object.

## References

Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about
Markov chains. \*The Annals of Mathematical Statistics\*, 28(1), 89–110.

Kullback, S., Kupperman, M. and Ku, H. H. (1962). Tests for Contingency
Tables and Markov Chains. \*Technometrics\*, 4(4), 573–608.

Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about
Markov chains. \*The Annals of Mathematical Statistics\*, 28(1), 89–110.
