# Function to fit a CTMC

This function fits the underlying CTMC give the state transition data
and the transition times using the maximum likelihood method (MLE)

## Usage

``` r
ctmcFit(data, byrow = TRUE, name = "", confidencelevel = 0.95)
```

## Arguments

- data:

  It is a list of two elements. The first element is a character vector
  denoting the states. The second is a numeric vector denoting the
  corresponding transition times.

- byrow:

  Determines if the output transition probabilities of the underlying
  embedded DTMC are by row.

- name:

  Optional name for the CTMC.

- confidencelevel:

  Confidence level for the confidence interval construnction.

## Value

It returns a list containing the CTMC object and the confidence
intervals.

## Details

Note that in data, there must exist an element wise corresponding
between the two elements of the list and that data\[\[2\]\]\[1\] is
always 0.

## References

Continuous Time Markov Chains (vignette), Sai Bhargav Yalamanchi,
Giorgio Alfredo Spedicato 2015

## See also

[`rctmc`](rctmc.md)

## Author

Sai Bhargav Yalamanchi

## Examples

``` r
data <- list(c("a", "b", "c", "a", "b", "a", "c", "b", "c"), c(0, 0.8, 2.1, 2.4, 4, 5, 5.9, 8.2, 9))
ctmcFit(data)
#> $estimate
#> An object of class "ctmc"
#> Slot "states":
#> [1] "a" "b" "c"
#> 
#> Slot "byrow":
#> [1] TRUE
#> 
#> Slot "generator":
#>            a          b          c
#> a -0.9090909  0.6060606  0.3030303
#> b  0.3225806 -0.9677419  0.6451613
#> c  0.3846154  0.3846154 -0.7692308
#> 
#> Slot "name":
#> [1] ""
#> 
#> 
#> $errors
#> $errors$dtmcConfidenceInterval
#> $errors$dtmcConfidenceInterval$confidenceLevel
#> [1] 0.95
#> 
#> $errors$dtmcConfidenceInterval$lowerEndpointMatrix
#>            a          b          c
#> a 0.00000000 0.20765960 0.06149194
#> b 0.06149194 0.00000000 0.20765960
#> c 0.09453121 0.09453121 0.00000000
#> 
#> $errors$dtmcConfidenceInterval$upperEndpointMatrix
#>           a         b         c
#> a 0.5614970 0.9385081 0.7923404
#> b 0.7923404 0.5614970 0.9385081
#> c 0.9054688 0.9054688 0.6576198
#> 
#> 
#> $errors$lambdaConfidenceInterval
#> $errors$lambdaConfidenceInterval$lowerEndpointVector
#> [1] 0.04576665 0.04871934 0.00000000
#> 
#> $errors$lambdaConfidenceInterval$upperEndpointVector
#> [1]  0.04576665  0.04871934 -0.12545166
#> 
#> 
#> 
```
