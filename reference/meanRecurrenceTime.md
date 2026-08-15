# Mean recurrence time

Computes the expected time to return to a recurrent state in case the
Markov chain starts there

## Usage

``` r
meanRecurrenceTime(object)
```

## Arguments

- object:

  the markovchain object

## Value

For a Markov chain it outputs is a named vector with the expected time
to first return to a state when the chain starts there. States present
in the vector are only the recurrent ones. If the matrix is ergodic
(i.e. irreducible), then all states are present in the output and order
is the same as states order for the Markov chain

## References

C. M. Grinstead and J. L. Snell. Introduction to Probability. American
Mathematical Soc., 2012.

## Author

Ignacio Cordón

## Examples

``` r
m <- matrix(1 / 10 * c(6,3,1,
                       2,3,5,
                       4,1,5), ncol = 3, byrow = TRUE)
mc <- new("markovchain", states = c("s","c","r"), transitionMatrix = m)
meanRecurrenceTime(mc)
#>        s        c        r 
#> 2.266667 4.250000 3.090909 
```
