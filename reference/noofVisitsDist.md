# Expected fraction of the first N steps spent in each state

Given the initial state \\i\\, returns for every state \\j\\ the
expected fraction of the first `N` steps that the DTMC spends in \\j\\.

## Usage

``` r
noofVisitsDist(markovchain,N,state)
```

## Arguments

- markovchain:

  a markovchain-class object

- N:

  number of steps, a positive integer

- state:

  the initial state

## Value

a named numeric vector with one element per state, summing to one.

## Details

The value for state \\j\\ is \$\$\frac{1}{N}\sum\_{k=1}^{N} (P^k)\_{ij}
= \frac{E\[V_j(N)\]}{N},\$\$ where \\V_j(N)\\ is the number of visits to
\\j\\ at times \\1, \dots, N\\ (the initial state, at time 0, is not
counted). The values sum to one, and multiplied by `N` they give the
expected numbers of visits. As `N` grows they converge to the stationary
distribution for an irreducible chain.

Despite the name of the function, and the title of earlier versions of
this page, the result is not the joint distribution of the numbers of
visits \\(V_1(N), \dots, V_n(N))\\, which the package does not compute
(see issue \#139).

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

# expected numbers of visits during the first 5 steps
5 * noofVisitsDist(simpleMc,5,"a")
#>       a       b 
#> 1.74074 3.25926 
```
