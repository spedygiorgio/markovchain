# Fit a higher order multivariate Markov chain

Fits the higher order multivariate Markov chain of Ching, Ng and Fung
(2008) to several categorical sequences observed together: the next
state of each sequence depends on the last `order` states of all the
sequences, through a weighted sum of lag-specific transition matrices.

## Usage

``` r
fitHighOrderMultivarMC(seqMat, order = 2, Norm = 2)
```

## Arguments

- seqMat:

  A matrix or data frame with one categorical sequence per column, all
  of the same length; numeric or factor values are taken as state labels
  (compared as character strings).

- order:

  Order of the chain (default 2).

- Norm:

  Norm used to estimate the weights (default 2): the weights minimize
  this norm of the difference between the state frequencies of each
  sequence and their image under the model, as in Ching et al. (2008).
  The minimization needs the suggested package Rsolnp.

## Value

An [`hommc`](hommc-class.md) object. If Rsolnp is not available, `NULL`
with a warning.

## References

Ching, W.-K., Ng, M. K. and Fung, E. S. (2008). Higher-order
multivariate Markov chains and their applications. Linear Algebra and
its Applications 428(2-3), 492-507.
[doi:10.1016/j.laa.2007.05.021](https://doi.org/10.1016/j.laa.2007.05.021)

## See also

[`predictHommc()`](predictHommc.md) to simulate from the fit.

Other higher order chains: [`fitHigherOrder()`](fitHigherOrder.md),
[`predictHommc()`](predictHommc.md)

## Author

Giorgio Spedicato, Deepak Yadav

## Examples

``` r
data <- matrix(c("2", "1", "3", "3", "4", "3", "2", "1", "3", "3", "2", "1",
                 "2", "4", "4", "4", "4", "2", "3", "3", "1", "4", "3", "3"),
               ncol = 2)
if (requireNamespace("Rsolnp", quietly = TRUE)) {
  fitHighOrderMultivarMC(data, order = 2, Norm = 2)
}
#> This function is experimental
#> Order of multivariate markov chain = 2 
#> states = 1 2 3 4 
#> 
#> List of Lambda's and the corresponding transition matrix (by cols) :
#> Lambda1(1,1) : 0.2496703
#> P1(1,1) : 
#>   1 2   3 4
#> 1 0 1 0.0 0
#> 2 0 0 0.4 0
#> 3 1 0 0.4 1
#> 4 0 0 0.2 0
#> 
#> Lambda2(1,1) : 6.559844e-05
#> P2(1,1) : 
#>   1 2   3 4
#> 1 0 0 0.4 0
#> 2 0 0 0.2 1
#> 3 1 1 0.2 0
#> 4 0 0 0.2 0
#> 
#> Lambda1(1,2) : 0.7501985
#> P1(1,2) : 
#>   1   2         3   4
#> 1 0 0.5 0.6666667 0.0
#> 2 0 0.5 0.0000000 0.2
#> 3 1 0.0 0.3333333 0.6
#> 4 0 0.0 0.0000000 0.2
#> 
#> Lambda2(1,2) : 6.559859e-05
#> P2(1,2) : 
#>   1   2 3   4
#> 1 0 0.5 0 0.2
#> 2 1 0.0 0 0.2
#> 3 0 0.5 1 0.4
#> 4 0 0.0 0 0.2
#> 
#> Lambda1(2,1) : 0.2692313
#> P1(2,1) : 
#>     1         2   3 4
#> 1 0.5 0.0000000 0.0 0
#> 2 0.0 0.0000000 0.0 1
#> 3 0.0 0.6666667 0.4 0
#> 4 0.5 0.3333333 0.6 0
#> 
#> Lambda2(2,1) : 0.2692313
#> P2(2,1) : 
#>   1   2   3 4
#> 1 0 0.5 0.0 0
#> 2 0 0.0 0.2 0
#> 3 0 0.0 0.6 1
#> 4 1 0.5 0.2 0
#> 
#> Lambda1(2,2) : 0.4615374
#> P1(2,2) : 
#>   1   2         3   4
#> 1 0 0.0 0.3333333 0.0
#> 2 0 0.0 0.0000000 0.2
#> 3 0 0.5 0.6666667 0.2
#> 4 1 0.5 0.0000000 0.6
#> 
#> Lambda2(2,2) : 8.013338e-09
#> P2(2,2) : 
#>   1   2   3   4
#> 1 0 0.0 0.5 0.0
#> 2 0 0.0 0.0 0.2
#> 3 1 0.5 0.0 0.4
#> 4 0 0.5 0.5 0.4
#> 
```
