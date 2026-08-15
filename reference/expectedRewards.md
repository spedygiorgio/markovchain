# Expected Rewards for a markovchain

Given a markovchain object and reward values for every state, function
calculates expected reward value after n steps.

## Usage

``` r
expectedRewards(markovchain,n,rewards)
```

## Arguments

- markovchain:

  the markovchain-class object

- n:

  no of steps of the process

- rewards:

  vector depicting rewards coressponding to states

## Value

returns a vector of expected rewards for different initial states

## Details

the function uses a dynamic programming approach to solve a recursive
equation described in reference.

## References

Stochastic Processes: Theory for Applications, Robert G. Gallager,
Cambridge University Press

## Author

Vandit Jain

## Examples

``` r
transMatr<-matrix(c(0.99,0.01,0.01,0.99),nrow=2,byrow=TRUE)
simpleMc<-new("markovchain", states=c("a","b"),
             transitionMatrix=transMatr)
expectedRewards(simpleMc,1,c(0,1))
#> [1] 0.01 1.99
```
