# Transition probability between two states

Returns the one-step probability of moving from state `t0` to state
`t1`.

## Usage

``` r
transitionProbability(object, t0, t1)

# S4 method for class 'markovchain'
transitionProbability(object, t0, t1)
```

## Arguments

- object:

  A [`markovchain`](markovchain-class.md) object.

- t0:

  Name of the current state.

- t1:

  Name of the next state.

## Value

A number, the transition probability.

## See also

[`conditionalDistribution()`](conditionalDistribution.md) for the whole
distribution of the next state.

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
transitionProbability(markovB,"b", "c")
#> [1] 0
```
