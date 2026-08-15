# A function to compute multinomial confidence intervals of DTMC

Return estimated transition matrix assuming a Multinomial Distribution

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

  An estimated transition matrix.

- countsTransitionMatrix:

  Empirical (conts) transition matrix, on which the `transitionMatrix`
  was performed.

- confidencelevel:

  confidence interval level.

## Value

Two matrices containing the confidence intervals.

## References

Constructing two-sided simultaneous confidence intervals for multinomial
proportions for small counts in a large number of cells. Journal of
Statistical Software 5(6) (2000)

## See also

`markovchainFit`

## Examples

``` r
seq<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
mcfit<-markovchainFit(data=seq,byrow=TRUE)
seqmat<-createSequenceMatrix(seq)
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
