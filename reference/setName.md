# Set the name of a Markov chain

Replaces the name of a [`markovchain`](markovchain-class.md) object (its
`name` slot).

## Usage

``` r
name(object) <- value

# S4 method for class 'markovchain'
name(object) <- value
```

## Arguments

- object:

  A [`markovchain`](markovchain-class.md) object.

- value:

  The new name, a character string.

## Value

The object with the new name (used as `name(object) <- value`).

## See also

[`name()`](getName.md) to read it.

## Author

Giorgio Spedicato, Deepak Yadav

## Examples

``` r
statesNames <- c("a", "b", "c")
markovB <- new("markovchain", states = statesNames, transitionMatrix =
                matrix(c(0.2, 0.5, 0.3, 0, 1, 0, 0.1, 0.8, 0.1), nrow = 3,
                byrow = TRUE, dimnames=list(statesNames,statesNames)),
                name = "A markovchain Object" 
)
name(markovB) <- "dangerous mc"
name(markovB)
#> [1] "dangerous mc"
```
