# Function to fit a discrete Markov chain

Given a sequence of states arising from a stationary state, it fits the
underlying Markov chain distribution using either MLE (also using a
Laplacian smoother), bootstrap or by MAP (Bayesian) inference.

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

  It can be a character vector or a \$\$n x n\$\$ matrix or a \$\$n x
  n\$\$ data frame or a list

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

  \$\$\alpha\$\$ level for conficence intervals width. Used only when
  `method` equal to "mle".

- confint:

  a boolean to decide whether to compute Confidence Interval or not.

- hyperparam:

  Hyperparameter matrix for the a priori distribution. If none is
  provided, default value of 1 is assigned to each parameter. This must
  be of size \$\$k x k\$\$ where k is the number of states in the chain
  and the values should typically be non-negative integers.

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

  Possible states which are not present in the given sequence

- progress:

  Should a text progress bar be shown? It is only used by the
  "bootstrap" method, the other methods being fast; see Details.

- stringchar:

  It can be a \$\$n x n\$\$ matrix or a character vector or a list

- toRowProbs:

  converts a sequence matrix into a probability matrix

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

A list containing an estimate, log-likelihood, and, when "bootstrap"
method is used, a matrix of standards deviations and the bootstrap
samples. When the "mle", "bootstrap" or "map" method is used, the lower
and upper confidence bounds are returned along with the standard error.
The "map" method also returns the expected value of the parameters with
respect to the posterior distribution.

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

A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall
2010

Inferring Markov Chains: Bayesian Estimation, Model Comparison, Entropy
Rate, and Out-of-Class Modeling, Christopher C. Strelioff, James P.
Crutchfield, Alfred Hubler, Santa Fe Institute

Yalamanchi SB, Spedicato GA (2015). Bayesian Inference of First Order
Markov Chains. R package version 0.2.5

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
