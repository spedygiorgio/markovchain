# Name of a Markov chain

Returns the name of a [`markovchain`](markovchain-class.md) object (its
`name` slot).

## Usage

``` r
name(object)

# S4 method for class 'markovchain'
name(object)
```

## Arguments

- object:

  A [`markovchain`](markovchain-class.md) object.

## Value

A character string.

## See also

`name<-` to change it.

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
name(markovB)
#> [1] "A markovchain Object"
```
