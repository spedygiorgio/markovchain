# Check whether a Markov chain is irreducible

A Markov chain is irreducible if all its states communicate, i.e. form a
single communicating class.

## Usage

``` r
is.irreducible(object)
```

## Arguments

- object:

  A [`markovchain`](markovchain-class.md) object.

## Value

`TRUE` if the chain is irreducible, `FALSE` otherwise.

## Details

Irreducibility does not imply aperiodicity: an irreducible chain on
finitely many states has a unique stationary distribution, but it
converges to it from every initial state only if it is also aperiodic
(see [`period()`](structuralAnalysis.md) and
[`is.regular()`](is.regular.md)).

## References

Feres, Matlab listings for Markov Chains.

## See also

[`communicatingClasses()`](structuralAnalysis.md),
[`is.regular()`](is.regular.md), [`period()`](structuralAnalysis.md)

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
