# `conditionalDistribution` of a Markov Chain

It extracts the conditional distribution of the subsequent state, given
current state.

## Usage

``` r
conditionalDistribution(object, state)
```

## Arguments

- object:

  A `markovchain` object.

- state:

  Subsequent state.

## Value

A named probability vector

## References

A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall
2010

## See also

[`markovchain`](markovchain-class.md)

## Author

Giorgio Spedicato, Deepak Yadav

## Examples

``` r
# define a markov chain
statesNames <- c("a", "b", "c")
markovB <- new("markovchain", states = statesNames, transitionMatrix = 
               matrix(c(0.2, 0.5, 0.3, 0, 1, 0, 0.1, 0.8, 0.1),nrow = 3, 
                      byrow = TRUE, dimnames = list(statesNames, statesNames)))
                      
conditionalDistribution(markovB, "b")                       
#> a b c 
#> 0 1 0 
```
