# Simulate from a homogeneous or non-homogeneous Markov chain

Simulates trajectories of a homogeneous Markov chain (a
[`markovchain`](markovchain-class.md) object) or of a non-homogeneous
one (a [`markovchainList`](markovchainList-class.md) object, whose
\\k\\-th matrix gives the transition from time \\k - 1\\ to time \\k\\).

## Usage

``` r
rmarkovchain(
  n,
  object,
  what = "data.frame",
  useRCpp = TRUE,
  parallel = FALSE,
  num.cores = NULL,
  ...
)
```

## Arguments

- n:

  For a `markovchain`, the number of steps of the trajectory; for a
  `markovchainList`, the number of trajectories, each as long as the
  list.

- object:

  A [`markovchain`](markovchain-class.md) or
  [`markovchainList`](markovchainList-class.md) object.

- what:

  Form of the output for a `markovchainList`: `"data.frame"` (default;
  columns `iteration` and `values`, one row per simulated state),
  `"matrix"` (one row per trajectory) or `"list"` (one character vector
  per trajectory). Ignored for a `markovchain`.

- useRCpp:

  Logical: use the C++ implementation (default)?

- parallel:

  Logical: simulate the trajectories of a `markovchainList` in parallel
  (default `FALSE`)?

- num.cores:

  Number of threads used when `parallel = TRUE`. If `NULL` (the default)
  the thread count is read from `getOption("RcppParallel.numThreads")`,
  then `getOption("Ncpus")`, then the environment variables
  `RCPP_PARALLEL_NUM_THREADS` and `OMP_NUM_THREADS`, falling back to
  `min(2, cores)` as CRAN policy requires.

- ...:

  Further arguments: `t0` (initial state) and `include.t0` (include it
  in the output), as in
  [`markovchainSequence()`](markovchainSequence.md).

## Value

For a `markovchain`, a character vector of states. For a
`markovchainList`, a data frame, matrix or list, according to `what`.

## Details

For a `markovchain` the function calls
[`markovchainSequence()`](markovchainSequence.md). For a
`markovchainList` each trajectory starts from `t0` (by default a state
drawn at random) and moves once with each matrix of the list.

## References

Ross, S. (2010). A First Course in Probability (8th edition). Prentice
Hall.

## See also

[`markovchainFit()`](markovchainFit.md) and
[`markovchainListFit()`](markovchainListFit.md) for the reverse
operation.

Other simulation: [`markovchainSequence()`](markovchainSequence.md)

## Author

Giorgio Spedicato

## Examples

``` r
statesNames <- c("a", "b", "c")
mcA <- new("markovchain", states = statesNames, transitionMatrix =
  matrix(c(0.2, 0.5, 0.3, 0, 0.2, 0.8, 0.1, 0.8, 0.1), nrow = 3,
  byrow = TRUE, dimnames = list(statesNames, statesNames)))
rmarkovchain(n = 20, object = mcA, t0 = "a")
#>  [1] "b" "b" "c" "b" "c" "b" "c" "b" "c" "b" "c" "b" "c" "b" "b" "c" "b" "b" "b"
#> [20] "c"

# a non-homogeneous chain: a different matrix at each of three steps
mcB <- new("markovchain", states = statesNames, transitionMatrix =
  matrix(c(0.5, 0.5, 0, 0.1, 0.8, 0.1, 0, 0.5, 0.5), nrow = 3,
  byrow = TRUE, dimnames = list(statesNames, statesNames)))
mclist <- new("markovchainList", markovchains = list(mcA, mcB, mcA))
rmarkovchain(5, mclist, what = "matrix")
#>      [,1] [,2] [,3]
#> [1,] "c"  "b"  "c" 
#> [2,] "b"  "b"  "b" 
#> [3,] "c"  "b"  "c" 
#> [4,] "b"  "b"  "c" 
#> [5,] "c"  "c"  "b" 
rmarkovchain(2, mclist, what = "list", t0 = "a", include.t0 = TRUE)
#> [[1]]
#> [1] "a" "c" "b" "c"
#> 
#> [[2]]
#> [1] "a" "b" "b" "b"
#> 
```
