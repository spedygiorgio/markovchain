# Function to get the transition probabilities from initial to subsequent states.

This is a convenience function to get transition probabilities.

## Usage

``` r
transitionProbability(object, t0, t1)

# S4 method for class 'markovchain'
transitionProbability(object, t0, t1)
```

## Arguments

- object:

  A `markovchain` object.

- t0:

  Initial state.

- t1:

  Subsequent state.

## Value

Numeric Vector

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
transitionProbability(markovB,"b", "c")
#> [1] 0
```
