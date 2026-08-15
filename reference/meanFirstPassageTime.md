# Mean First Passage Time for irreducible Markov chains

Given an irreducible (ergodic) markovchain object, this function
calculates the expected number of steps to reach other states

## Usage

``` r
meanFirstPassageTime(object, destination)
```

## Arguments

- object:

  the markovchain object

- destination:

  a character vector representing the states respect to which we want to
  compute the mean first passage time. Empty by default

## Value

a Matrix of the same size with the average first passage times if
destination is empty, a vector if destination is not

## Details

For an ergodic Markov chain it computes:

- If destination is empty, the average first time (in steps) that takes
  the Markov chain to go from initial state i to j. (i, j) represents
  that value in case the Markov chain is given row-wise, (j, i) in case
  it is given col-wise.

- If destination is not empty, the average time it takes us from the
  remaining states to reach the states in `destination`

## References

C. M. Grinstead and J. L. Snell. Introduction to Probability. American
Mathematical Soc., 2012.

## Author

Toni Giorgino, Ignacio Cordón

## Examples

``` r
m <- matrix(1 / 10 * c(6,3,1,
                       2,3,5,
                       4,1,5), ncol = 3, byrow = TRUE)
mc <- new("markovchain", states = c("s","c","r"), transitionMatrix = m)
meanFirstPassageTime(mc, "r")
#>        s        c 
#> 4.545455 2.727273 


# Grinstead and Snell's "Oz weather" worked out example
mOz <- matrix(c(2,1,1,
                2,0,2,
                1,1,2)/4, ncol = 3, byrow = TRUE)

mcOz <- new("markovchain", states = c("s", "c", "r"), transitionMatrix = mOz)
meanFirstPassageTime(mcOz)
#>          s c        r
#> s 0.000000 4 3.333333
#> c 2.666667 0 2.666667
#> r 3.333333 4 0.000000
```
