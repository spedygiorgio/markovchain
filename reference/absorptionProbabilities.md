# Absorption probabilities

Computes the absorption probability from each transient state to each
recurrent one (i.e. the (i, j) entry or (j, i), in a stochastic matrix
by columns, represents the probability that the first not transient
state we can go from the transient state i is j (and therefore we are
going to be absorbed in the communicating recurrent class of j)

## Usage

``` r
absorptionProbabilities(object)
```

## Arguments

- object:

  the markovchain object

## Value

A named vector with the expected number of steps to go from a transient
state to any of the recurrent ones

## References

C. M. Grinstead and J. L. Snell. Introduction to Probability. American
Mathematical Soc., 2012.

## Author

Ignacio Cordón

## Examples

``` r
m <- matrix(c(1/2, 1/2, 0,
              1/2, 1/2, 0,
                0, 1/2, 1/2), ncol = 3, byrow = TRUE)
mc <- new("markovchain", states = letters[1:3], transitionMatrix = m)
absorptionProbabilities(mc)
#>   a b
#> c 0 1
```
