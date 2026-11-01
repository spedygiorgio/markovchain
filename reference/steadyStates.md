# Stationary distributions of a Markov chain

Returns the stationary distributions of a Markov chain: the probability
vectors \\\pi\\ with \\\pi P = \pi\\. There is one for each recurrent
class; every stationary distribution is a convex combination of them.

## Usage

``` r
steadyStates(object)
```

## Arguments

- object:

  A [`markovchain`](markovchain-class.md) object.

## Value

A matrix with one stationary distribution per row (per column if the
chain is stored by columns, `byrow = FALSE`), with the states as column
(row) names.

## Note

The steady states are identified starting from which eigenvectors
correspond to identity eigenvalues and then normalizing them to sum up
to unity. When negative values are found in the matrix, the eigenvalues
extraction is performed on the recurrent classes submatrix.

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
steadyStates(markovB)
#>      a b c
#> [1,] 0 1 0
```
