# Functions to fit a higher order Markov chain

Given a sequence of states arising from a stationary state, it fits the
underlying Markov chain distribution with higher order.

## Usage

``` r
fitHigherOrder(sequence, order = 2)
seq2freqProb(sequence)
seq2matHigh(sequence, order)
```

## Arguments

- sequence:

  A character list.

- order:

  Markov chain order

## Value

A list containing lambda, Q, and X.

## References

Ching, W. K., Huang, X., Ng, M. K., & Siu, T. K. (2013). Higher-order
markov chains. In Markov Chains (pp. 141-176). Springer US.

Ching, W. K., Ng, M. K., & Fung, E. S. (2008). Higher-order multivariate
Markov chains and their applications. Linear Algebra and its
Applications, 428(2), 492-507.

## Author

Giorgio Spedicato, Tae Seung Kang

## Examples

``` r
sequence<-c("a", "a", "b", "b", "a", "c", "b", "a", "b", "c", "a", "b",
            "c", "a", "b", "c", "a", "b", "a", "b")
fitHigherOrder(sequence)
#> $lambda
#> [1] 1.000000e+00 1.626303e-08
#> 
#> $Q
#> $Q[[1]]
#>       a         b    c
#> a 0.125 0.4285714 0.75
#> b 0.750 0.1428571 0.25
#> c 0.125 0.4285714 0.00
#> 
#> $Q[[2]]
#>           a         b    c
#> a 0.1428571 0.5714286 0.25
#> b 0.4285714 0.2857143 0.75
#> c 0.4285714 0.1428571 0.00
#> 
#> 
#> $X
#>   a   b   c 
#> 0.4 0.4 0.2 
#> 
```
