# Fit a discrete time Markov chain

Estimates the transition matrix of a homogeneous discrete time Markov
chain from observed sequences of states, by maximum likelihood
(optionally with Laplace smoothing), by bootstrap, or by maximum a
posteriori (Bayesian) estimation with Dirichlet priors.
`createSequenceMatrix()` returns the matrix of transition counts (or
frequencies) on which the estimates are based.

## Usage

``` r
.markovchainFitRcpp(
  data,
  method = "mle",
  byrow = TRUE,
  nboot = 10L,
  laplacian = 0,
  name = "",
  parallel = FALSE,
  confidencelevel = 0.95,
  confint = TRUE,
  hyperparam = matrix(),
  sanitize = FALSE,
  possibleStates = character(),
  progress = FALSE
)

createSequenceMatrix(
  stringchar,
  toRowProbs = FALSE,
  sanitize = FALSE,
  possibleStates = character()
)

markovchainFit(data, method = "mle", byrow = TRUE, nboot = 10L,
  laplacian = 0, name = "", parallel = FALSE, confidencelevel = 0.95,
  confint = TRUE, hyperparam = matrix(), sanitize = FALSE,
  possibleStates = character(), absorbingStates = character(),
  progress = FALSE, num.cores = NULL)
```

## Arguments

- data:

  The observed sequence(s): a character vector (one sequence), a list of
  character vectors (independent sequences), or a matrix or data frame
  with one trajectory per row (per column if `byrow = FALSE`) and one
  time per column.

- method:

  Method used to estimate the Markov chain. Either "mle", "map",
  "bootstrap" or "laplace". All four are available for a single
  sequence. For a list of sequences, "mle", "map" and "laplace" pool the
  transition counts over the sequences, while "bootstrap" is not
  available and raises an error explaining why.

- byrow:

  For a character vector or a list of sequences, it tells whether the
  fitted transition matrix is stored by row (the default) or by column;
  the `byrow` slot of the returned chain records the choice. For matrix
  or data frame input it instead describes the input data – whether each
  observed trajectory is a row or a column of `data` – and the fitted
  chain is stored by row either way.

- nboot:

  Number of bootstrap replicates in case "bootstrap" is used.

- laplacian:

  Laplacian smoothing parameter, default zero. It is only used when
  "laplace" method is chosen.

- name:

  Optional character for name slot.

- parallel:

  Use parallel processing when performing Boostrap estimates.

- confidencelevel:

  Confidence level of the intervals (default 0.95).

- confint:

  Logical: compute standard errors and confidence intervals (default
  `TRUE`)?

- hyperparam:

  For `method = "map"`, the matrix of Dirichlet hyperparameters, \\k
  \times k\\ for \\k\\ states, with the states as dimnames (see
  [`inferHyperparam()`](inferHyperparam.md)). The default sets all of
  them to 1.

- sanitize:

  how to deal with the states that have no observed outgoing transition,
  which is what `possibleStates` typically introduces. `FALSE` (the
  default) leaves their row at zero, which makes a transition matrix
  that is not stochastic. `TRUE`, or equivalently `"uniform"`, puts 1 in
  every entry of such a row, so that the row becomes a uniform
  distribution over all the states: the unobserved states are then
  assumed to move to any state with equal probability, which is an
  assumption about the data and not a consequence of it. `"absorbing"`
  instead puts 1 on the diagonal only, making every unobserved state
  absorbing; this also gives a stochastic matrix, but adds no transition
  that was never observed (see \#213).

- possibleStates:

  States to include in the result even if they do not appear in the
  data.

- progress:

  Should a text progress bar be shown? It is only used by the
  "bootstrap" method, the other methods being fast; see Details.

- stringchar:

  For `createSequenceMatrix()`: a character vector (one sequence), a
  list of sequences, or a two-column matrix of transitions (from, to),
  one per row.

- toRowProbs:

  Logical: return the transition frequencies (each row divided by its
  total) instead of the counts?

- absorbingStates:

  Character vector of states that are known a priori to be absorbing.
  The corresponding rows are set to the identity row after MLE fitting
  when `byrow = TRUE`; the corresponding columns are set to the identity
  column when `byrow = FALSE`. The argument is currently supported only
  for `method = "mle"`.

- num.cores:

  Number of threads the parallel bootstrap path uses when
  `method = "bootstrap"` and `parallel = TRUE`. If `NULL` (the default)
  the thread count is read from `getOption("RcppParallel.numThreads")` /
  `getOption("Ncpus")` / `OMP_NUM_THREADS` /
  `RCPP_PARALLEL_NUM_THREADS`, falling back to `min(2, cores)` as CRAN
  policy requires. Ignored when `parallel = FALSE`.

## Value

`markovchainFit()` returns a list whose element `estimate` is the fitted
[`markovchain`](markovchain-class.md) object. The other elements depend
on the method:

- `"mle"`: `standardError`, `confidenceLevel`, `lowerEndpointMatrix`,
  `upperEndpointMatrix` (when `confint = TRUE`) and `logLikelihood`;

- `"laplace"`: `logLikelihood`;

- `"map"`: `expectedValue` (posterior mean of the transition matrix),
  `standardError`, `confidenceInterval` and `logLikelihood`;

- `"bootstrap"`: `standardError`, `confidenceInterval`,
  `bootStrapSamples` and `logLikelihood`.

For matrix or data frame input only the maximum likelihood estimate is
computed, whatever `method`, and there is no `logLikelihood`.

`createSequenceMatrix()` returns the square matrix of transition counts
(frequencies, with `toRowProbs = TRUE`), with the states as dimnames.

## Details

Disabling confint would lower the computation time on large datasets. If
`data` or `stringchar` contain `NAs`, the related `NA` containing
transitions will be ignored.

With `progress = TRUE` the "bootstrap" method shows a text progress bar
([`txtProgressBar`](https://rdrr.io/r/utils/txtProgressBar.html), style
3) covering the simulation of the `nboot` bootstrap sequences and the
estimation of a transition matrix from each of them. With
`parallel = TRUE` the sequences are simulated in parallel threads, which
cannot report progress, so the bar only covers the estimation step. The
computation can be interrupted while the bar is shown.

When `absorbingStates` is supplied, the declared states must have no
observed outgoing transitions. This allows terminal states in censored
customer journeys to be represented as absorbing states without adding
artificial observations.

`sanitize = "absorbing"` reaches the same result without naming the
states: every state that has no observed outgoing transition, which is
what `possibleStates` typically introduces, is made absorbing. Unlike
`absorbingStates`, it works with every `method`, since it only replaces
the uniform row that `sanitize = TRUE` would have produced. As with
`absorbingStates`, only the estimate is constrained: any confidence
bounds and standard errors keep the values the unconstrained fit
assigned to those rows.

## Note

This function has been rewritten in Rcpp. Bootstrap algorithm has been
defined "heuristically". In addition, parallel facility is not complete,
involving only a part of the bootstrap process. When `data` is either a
`data.frame` or a `matrix` object, only MLE fit is currently available.

## References

Anderson, T. W. and Goodman, L. A. (1957). Statistical inference about
Markov chains. The Annals of Mathematical Statistics 28(1), 89-110.

Strelioff, C. C., Crutchfield, J. P. and Hubler, A. W. (2007). Inferring
Markov chains: Bayesian estimation, model comparison, entropy rate, and
out-of-class modeling. Physical Review E 76, 011106.
[doi:10.1103/PhysRevE.76.011106](https://doi.org/10.1103/PhysRevE.76.011106)

Yalamanchi, S. B. and Spedicato, G. A. (2015). Bayesian inference of
first order Markov chains. Vignette of the markovchain package.

## See also

[`markovchainSequence`](markovchainSequence.md),
[`markovchainListFit`](markovchainListFit.md)

## Author

Giorgio Spedicato, Tae Seung Kang, Sai Bhargav Yalamanchi

## Examples

``` r
sequence <- c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", 
              "b", "b", "b", "a")        
sequenceMatr <- createSequenceMatrix(sequence, sanitize = FALSE)
mcFitMLE <- markovchainFit(data = sequence)
mcFitBSP <- markovchainFit(data = sequence, method = "bootstrap", nboot = 5, name = "Bootstrap Mc")

na.sequence <- c("a", NA, "a", "b")
# There will be only a (a,b) transition        
na.sequenceMatr <- createSequenceMatrix(na.sequence, sanitize = FALSE)
mcFitMLE <- markovchainFit(data = na.sequence)

# data can be a list of character vectors
sequences <- list(x = c("a", "b", "a"), y = c("b", "a", "b", "a", "c"))
mcFitMap <- markovchainFit(sequences, method = "map")
mcFitMle <- markovchainFit(sequences, method = "mle")
```
