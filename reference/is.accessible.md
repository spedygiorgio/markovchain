# Check whether a state can be reached from another

Checks whether state `to` is accessible from state `from`, i.e. whether
the chain started in `from` reaches `to` with positive probability in
some number of steps.

## Usage

``` r
is.accessible(object, from, to)
```

## Arguments

- object:

  A [`markovchain`](markovchain-class.md) object.

- from:

  Name of the starting state.

- to:

  Name of the target state.

## Value

`TRUE` if `to` is accessible from `from`, `FALSE` otherwise.

## References

James Montgomery, University of Madison

## See also

[`is.irreducible()`](is.irreducible.md),
[`communicatingClasses()`](structuralAnalysis.md)

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
