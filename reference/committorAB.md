# Calculates committor of a markovchain object with respect to set A, B

Returns the probability of hitting states rom set A before set B with
different initial states

## Usage

``` r
committorAB(object,A,B,p)
```

## Arguments

- object:

  a markovchain class object

- A:

  a set of states

- B:

  a set of states

- p:

  initial state (default value : 1)

## Value

Return a vector of probabilities in case initial state is not provided
else returns a number

## Details

The function solves a system of linear equations to calculate probaility
that the process hits a state from set A before any state from set B

## Examples

``` r
transMatr <- matrix(c(0,0,0,1,0.5,
                      0.5,0,0,0,0,
                      0.5,0,0,0,0,
                      0,0.2,0.4,0,0,
                      0,0.8,0.6,0,0.5),
                      nrow = 5)
object <- new("markovchain", states=c("a","b","c","d","e"),transitionMatrix=transMatr)
committorAB(object,c(5),c(3))
#>         a         b         c         d         e 
#> 0.4444444 0.8888889 0.0000000 0.4444444 1.0000000 
```
