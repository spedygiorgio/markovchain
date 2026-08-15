# Mean absorption time

Computes the expected number of steps to go from any of the transient
states to any of the recurrent states. The Markov chain should have at
least one transient state for this method to work

## Usage

``` r
meanAbsorptionTime(object)
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
times <- meanAbsorptionTime(mc)
```
