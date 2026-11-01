# Simulate a trajectory of a Markov chain

Simulates `n` steps of a homogeneous discrete time Markov chain from the
initial state `t0`, drawing each state from the transition probabilities
of the previous one.

## Usage

``` r
markovchainSequence(
  n,
  markovchain,
  t0 = sample(markovchain@states, 1),
  include.t0 = FALSE,
  useRCpp = TRUE
)
```

## Arguments

- n:

  Number of steps, a non-negative integer.

- markovchain:

  A [`markovchain`](markovchain-class.md) object.

- t0:

  Name of the initial state. The default draws it uniformly at random
  among the states.

- include.t0:

  Logical: include the initial state in the output? If `TRUE` the output
  has `n + 1` states.

- useRCpp:

  Logical: use the C++ implementation (default) or the R one? Both give
  the same sequence for the same random seed.

## Value

A character vector with the simulated states.

## References

Ross, S. (2010). A First Course in Probability (8th edition). Prentice
Hall.

## See also

[`rmarkovchain()`](rmarkovchain.md), which also simulates from
non-homogeneous chains ([`markovchainList`](markovchainList-class.md)
objects); [`markovchainFit()`](markovchainFit.md) for the reverse
operation.

Other simulation: [`rmarkovchain()`](rmarkovchain.md)

## Author

Giorgio Spedicato

## Examples

``` r
statesNames <- c("a", "b", "c")
mcB <- new("markovchain", states = statesNames,
  transitionMatrix = matrix(c(0.2, 0.5, 0.3, 0, 0.2, 0.8, 0.1, 0.8, 0.1),
  nrow = 3, byrow = TRUE, dimnames = list(statesNames, statesNames)))
markovchainSequence(n = 20, markovchain = mcB, t0 = "a")
#>  [1] "b" "c" "b" "c" "b" "c" "b" "c" "b" "b" "c" "b" "c" "c" "b" "c" "b" "c" "b"
#> [20] "c"
markovchainSequence(n = 20, markovchain = mcB, t0 = "a", include.t0 = TRUE)
#>  [1] "a" "b" "c" "b" "c" "a" "c" "b" "c" "b" "c" "b" "c" "b" "c" "b" "c" "a" "b"
#> [20] "c" "a"
```
