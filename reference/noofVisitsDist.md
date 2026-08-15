# return a joint pdf of the number of visits to the various states of the DTMC

This function would return a joint pdf of the number of visits to the
various states of the DTMC during the first N steps.

## Usage

``` r
noofVisitsDist(markovchain,N,state)
```

## Arguments

- markovchain:

  a markovchain-class object

- N:

  no of steps

- state:

  the initial state

## Value

a numeric vector depicting the above described probability density
function.

## Details

This function would return a joint pdf of the number of visits to the
various states of the DTMC during the first N steps.

## Author

Vandit Jain

## Examples

``` r
transMatr<-matrix(c(0.4,0.6,.3,.7),nrow=2,byrow=TRUE)
simpleMc<-new("markovchain", states=c("a","b"),
             transitionMatrix=transMatr, 
             name="simpleMc")   
noofVisitsDist(simpleMc,5,"a")
#>        a        b 
#> 0.348148 0.651852 
```
