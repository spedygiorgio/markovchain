# Hitting probabilities for markovchain

Given a markovchain object, this function calculates the probability of
ever arriving from state i to j

## Usage

``` r
hittingProbabilities(object)
```

## Arguments

- object:

  the markovchain-class object

## Value

a matrix of hitting probabilities

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
```
