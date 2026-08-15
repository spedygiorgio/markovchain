# function to calculate first passage probabilities

The function calculates first passage probability for a subset of states
given an initial state.

## Usage

``` r
firstPassageMultiple(object, state, set, n)
```

## Arguments

- object:

  a markovchain-class object

- state:

  intital state of the process (charactervector)

- set:

  set of states A, first passage of which is to be calculated

- n:

  Number of rows on which compute the distribution

## Value

A vector of size n showing the first time proabilities

## References

Renaldo Feres, Notes for Math 450 Matlab listings for Markov chains; MIT
OCW, course - 6.262, Discrete Stochastic Processes, course-notes, chap
-05

## See also

[`firstPassage`](firstPassage.md)

## Author

Vandit Jain

## Examples

``` r
statesNames <- c("a", "b", "c")
markovB <- new("markovchain", states = statesNames, transitionMatrix =
matrix(c(0.2, 0.5, 0.3,
         0, 1, 0,
         0.1, 0.8, 0.1), nrow = 3, byrow = TRUE,
       dimnames = list(statesNames, statesNames)
))
firstPassageMultiple(markovB,"a",c("b","c"),4)  
#>      set
#> 1 0.8000
#> 2 0.4000
#> 3 0.1190
#> 4 0.0379
```
