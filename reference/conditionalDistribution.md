# Conditional distribution of the next state

Returns the distribution of the next state of a Markov chain given its
current state: the row of the transition matrix of that state (its
column if the chain is stored by columns).

## Usage

``` r
conditionalDistribution(object, state)
```

## Arguments

- object:

  A [`markovchain`](markovchain-class.md) object.

- state:

  Name of the current state.

## Value

A probability vector named by the states.

## See also

[`transitionProbability()`](transitionProbability.md) for a single
probability.

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
