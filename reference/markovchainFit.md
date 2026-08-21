# Function to fit a discrete Markov chain

Given a sequence of states arising from a stationary state, it fits the
underlying Markov chain distribution using either MLE (also using a
Laplacian smoother), bootstrap or by MAP (Bayesian) inference.

## Usage

``` r
createSequenceMatrix(
  stringchar,
  toRowProbs = FALSE,
  sanitize = FALSE,
  possibleStates = character()
)

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
  possibleStates = character()
)

markovchainFit(data, method = "mle", byrow = TRUE, nboot = 10L,
  laplacian = 0, name = "", parallel = FALSE, confidencelevel = 0.95,
  confint = TRUE, hyperparam = matrix(), sanitize = FALSE,
  possibleStates = character(), absorbingStates = character())
```

## Arguments

- stringchar:

  It can be a \$\$n x n\$\$ matrix or a character vector or a list

- toRowProbs:

  converts a sequence matrix into a probability matrix

- sanitize:

  put 1 in all rows having rowSum equal to zero

- possibleStates:

  Possible states which are not present in the given sequence

- data:

  It can be a character vector or a \$\$n x n\$\$ matrix or a \$\$n x
  n\$\$ data frame or a list

- method:

  Method used to estimate the Markov chain. Either "mle", "map",
  "bootstrap" or "laplace"

- byrow:

  it tells whether the output Markov chain should show the transition
  probabilities by row.

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

- absorbingStates:

  Character vector of states that are known a priori to be absorbing.
  The corresponding rows are set to the identity row after MLE fitting
  when `byrow = TRUE`; the corresponding columns are set to the identity
  column when `byrow = FALSE`. The argument is currently supported only
  for `method = "mle"`.

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

When `absorbingStates` is supplied, the declared states must have no
observed outgoing transitions. This allows terminal states in censored
customer journeys to be represented as absorbing states without adding
artificial observations.

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
