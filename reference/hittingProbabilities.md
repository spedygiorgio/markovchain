# Hitting probabilities for markovchain

Given a markovchain object, this function calculates the probability of
ever arriving from state i to j

## Usage

``` r
hittingProbabilities(object, targets = NULL)
```

## Arguments

- object:

  the markovchain-class object

- targets:

  optional character vector of state names: only the hitting
  probabilities *towards* these states are computed. The default,
  `NULL`, means all the states, which gives the full matrix as before.
  Every target is handled independently of the others, so the work is
  proportional to the number of targets: for a large chain, asking only
  for the states of interest is much faster than computing the whole
  matrix and subsetting it. Duplicated or unknown names are an error.

## Value

a matrix of hitting probabilities. Entry `[i, j]` is the probability of
ever arriving from state `i` to state `j` (the probability of returning,
after at least one transition, on the diagonal); for a chain with
`byrow = FALSE` the matrix is transposed, as the transition matrix is.
With `targets`, only the columns (rows if `byrow = FALSE`) of the
targets are returned, in the order given, and they coincide with those
of the full matrix.

## References

R. Vélez, T. Prieto, Procesos Estocásticos, Librería UNED, 2013

## Author

Ignacio Cordón

## Examples

``` r
M <- markovchain:::zeros(5)
M[1,1] <- M[5,5] <- 1
M[2,1] <- M[2,3] <- 1/2
M[3,2] <- M[3,4] <- 1/2
M[4,2] <- M[4,5] <- 1/2

mc <- new("markovchain", transitionMatrix = M)
hittingProbabilities(mc)
#>     1     2     3         4   5
#> 1 1.0 0.000 0.000 0.0000000 0.0
#> 2 0.8 0.375 0.500 0.3333333 0.2
#> 3 0.6 0.750 0.375 0.6666667 0.4
#> 4 0.4 0.500 0.250 0.1666667 0.6
#> 5 0.0 0.000 0.000 0.0000000 1.0

# only the probabilities of ever reaching the first state
hittingProbabilities(mc, targets = "1")
#>     1
#> 1 1.0
#> 2 0.8
#> 3 0.6
#> 4 0.4
#> 5 0.0
```
