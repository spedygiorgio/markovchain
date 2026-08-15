# Defined states of a transition matrix

This method returns the states of a transition matrix.

## Usage

``` r
states(object)

# S4 method for class 'markovchain'
states(object)
```

## Arguments

- object:

  A discrete `markovchain` object

## Value

The character vector corresponding to states slot.

## References

A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall
2010

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
