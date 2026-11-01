# States of a Markov chain

Returns the names of the states of a Markov chain, in the order of its
transition matrix.

## Usage

``` r
states(object)

# S4 method for class 'markovchain'
states(object)
```

## Arguments

- object:

  A [`markovchain`](markovchain-class.md) object (or a
  [`ctmc`](ctmc-class.md) object).

## Value

A character vector, the `states` slot.

## See also

[`markovchain`](markovchain-class.md)

## Author

Giorgio Spedicato

## Examples

``` r
statesNames <- c("a", "b", "c")
markovB <- new("markovchain", states = statesNames, transitionMatrix =
                matrix(c(0.2, 0.5, 0.3, 0, 1, 0, 0.1, 0.8, 0.1), nrow = 3,
                byrow = TRUE, dimnames=list(statesNames,statesNames)),
                name = "A markovchain Object" 
)
states(markovB)
#> [1] "a" "b" "c"
names(markovB)
#> [1] "a" "b" "c"
```
