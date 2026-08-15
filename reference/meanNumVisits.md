# Mean num of visits for markovchain, starting at each state

Given a markovchain object, this function calculates a matrix where the
element (i, j) represents the expect number of visits to the state j if
the chain starts at i (in a Markov chain by columns it would be the
element (j, i) instead)

## Usage

``` r
meanNumVisits(object)
```

## Arguments

- object:

  the markovchain-class object

## Value

a matrix with the expect number of visits to each state

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
meanNumVisits(mc)
#>     1   2   3   4   5
#> 1 Inf 0.0 0.0 0.0   0
#> 2 Inf 0.6 0.8 0.4 Inf
#> 3 Inf 1.2 0.6 0.8 Inf
#> 4 Inf 0.8 0.4 0.2 Inf
#> 5   0 0.0 0.0 0.0 Inf
```
