# Return the generator matrix for a corresponding transition matrix

Calculate the generator matrix for a corresponding transition matrix

## Usage

``` r
transition2Generator(P, t = 1, method = "logarithm")
```

## Arguments

- P:

  transition matrix between time 0 and t

- t:

  time of observation

- method:

  "logarithm" returns the Matrix logarithm of the transition matrix

## Value

A matrix that represent the generator of P

## See also

[`rctmc`](rctmc.md)

## Examples

``` r
mymatr <- matrix(c(.4, .6, .1, .9), nrow = 2, byrow = TRUE)
Q <- transition2Generator(P = mymatr)
expm::expm(Q)
#>      [,1] [,2]
#> [1,]  0.4  0.6
#> [2,]  0.1  0.9
 
```
