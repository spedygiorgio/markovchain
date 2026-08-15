# Function to check if a Markov chain is irreducible (i.e. ergodic)

This function verifies whether a `markovchain` object transition matrix
is composed by only one communicating class.

## Usage

``` r
is.irreducible(object)
```

## Arguments

- object:

  A `markovchain` object

## Value

A boolean values.

## Details

It is based on `.communicatingClasses` internal function.

## References

Feres, Matlab listings for Markov Chains.

## See also

[`summary`](https://rdrr.io/r/base/summary.html)

## Author

Giorgio Spedicato

## Examples

``` r
statesNames <- c("a", "b")
mcA <- new("markovchain", transitionMatrix = matrix(c(0.7,0.3,0.1,0.9),
                                             byrow = TRUE, nrow = 2, 
                                             dimnames = list(statesNames, statesNames)
           ))
is.irreducible(mcA)
#> [1] TRUE
```
