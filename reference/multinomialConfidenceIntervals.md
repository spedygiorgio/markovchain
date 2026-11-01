# Simultaneous confidence intervals for transition probabilities

Simultaneous confidence intervals for the rows of an estimated
transition matrix, each row being a multinomial proportion, by the
method of Sison and Glaz (1995) in the implementation of May and Johnson
(2000).

## Usage

``` r
multinomialConfidenceIntervals(
  transitionMatrix,
  countsTransitionMatrix,
  confidencelevel = 0.95
)
```

## Arguments

- transitionMatrix:

  The estimated transition matrix.

- countsTransitionMatrix:

  The matrix of transition counts from which `transitionMatrix` was
  estimated (see [`createSequenceMatrix()`](markovchainFit.md)).

- confidencelevel:

  Confidence level (default 0.95).

## Value

A list with `confidenceLevel`, `lowerEndpointMatrix` and
`upperEndpointMatrix`, the latter two with the dimnames of
`transitionMatrix`.

## References

Sison, C. P. and Glaz, J. (1995). Simultaneous confidence intervals and
sample size determination for multinomial proportions. Journal of the
American Statistical Association 90, 366-369.

May, W. L. and Johnson, W. D. (2000). Constructing two-sided
simultaneous confidence intervals for multinomial proportions for small
counts in a large number of cells. Journal of Statistical Software 5(6).

## See also

[`markovchainFit()`](markovchainFit.md), whose `method = "mle"` gives
intervals for each probability separately rather than simultaneous ones
for each row.

## Examples

``` r
seq <- c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a",
         "b", "b", "b", "a")
mcfit <- markovchainFit(data = seq, byrow = TRUE)
seqmat <- createSequenceMatrix(seq)
multinomialConfidenceIntervals(mcfit$estimate@transitionMatrix, seqmat, 0.95)
#> $confidenceLevel
#> [1] 0.95
#> 
#> $lowerEndpointMatrix
#>           a         b
#> a 0.2222222 0.3333333
#> b 0.5714286 0.1428571
#> 
#> $upperEndpointMatrix
#>           a         b
#> a 0.8111456 0.9222567
#> b 1.0000000 0.6839473
#> 
```
