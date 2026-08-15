# Function to generate a sequence of states from homogeneous Markov chains.

Provided any `markovchain` object, it returns a sequence of states
coming from the underlying stationary distribution.

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

  Sample size

- markovchain:

  `markovchain` object

- t0:

  The initial state

- include.t0:

  Specify if the initial state shall be used

- useRCpp:

  Boolean. Should RCpp fast implementation being used? Default is yes.

## Value

A Character Vector

## Details

A sequence of size n is sampled.

## References

A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall
2010

## See also

[`markovchainFit`](markovchainFit.md)

## Author

Giorgio Spedicato

## Examples

``` r
# define the markovchain object
statesNames <- c("a", "b", "c")
mcB <- new("markovchain", states = statesNames, 
   transitionMatrix = matrix(c(0.2, 0.5, 0.3, 0, 0.2, 0.8, 0.1, 0.8, 0.1), 
   nrow = 3, byrow = TRUE, dimnames = list(statesNames, statesNames)))

# show the sequence
outs <- markovchainSequence(n = 100, markovchain = mcB, t0 = "a")
```
