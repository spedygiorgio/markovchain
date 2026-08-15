# Calculating full conditional probability using lower rate transition matrix

This function calculates full conditional probability at given time s
using lower rate transition matrix

## Usage

``` r
impreciseProbabilityatT(C,i,t,s,error,useRCpp)
```

## Arguments

- C:

  a ictmc class object

- i:

  initial state at time t

- t:

  initial time t. Default value = 0

- s:

  final time

- error:

  error rate. Default value = 0.001

- useRCpp:

  logical whether to use RCpp implementation; by default TRUE

## References

Imprecise Continuous-Time Markov Chains, Thomas Krak et al., 2016

## Author

Vandit Jain

## Examples

``` r
states <- c("n","y")
Q <- matrix(c(-1,1,1,-1),nrow = 2,byrow = TRUE,dimnames = list(states,states))
range <- matrix(c(1/52,3/52,1/2,2),nrow = 2,byrow = 2)
name <- "testictmc"
ictmc <- new("ictmc",states = states,Q = Q,range = range,name = name)
impreciseProbabilityatT(ictmc,2,0,1,10^-3,TRUE)
#> [1] 0.008259774 0.140983476
```
