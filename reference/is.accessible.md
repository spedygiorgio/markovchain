# Verify if a state j is reachable from state i.

This function verifies if a state is reachable from another, i.e., if
there exists a path that leads to state j leaving from state i with
positive probability

## Usage

``` r
is.accessible(object, from, to)
```

## Arguments

- object:

  A `markovchain` object.

- from:

  The name of state "i" (beginning state).

- to:

  The name of state "j" (ending state).

## Value

A boolean value.

## Details

It wraps an internal function named `reachabilityMatrix`.

## References

James Montgomery, University of Madison

## See also

`is.irreducible`

## Author

Giorgio Spedicato, Ignacio Cordón

## Examples

``` r
statesNames <- c("a", "b", "c")
markovB <- new("markovchain", states = statesNames, 
               transitionMatrix = matrix(c(0.2, 0.5, 0.3,
                                             0,   1,   0,
                                           0.1, 0.8, 0.1), nrow = 3, byrow = TRUE, 
                                         dimnames = list(statesNames, statesNames)
                                        )
               )
is.accessible(markovB, "a", "c")
#> [1] TRUE
```
