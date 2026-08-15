# Function to fit Higher Order Multivariate Markov chain

Given a matrix of categorical sequences it fits Higher Order
Multivariate Markov chain.

## Usage

``` r
fitHighOrderMultivarMC(seqMat, order = 2, Norm = 2)
```

## Arguments

- seqMat:

  a matrix or a data frame where each column is a categorical sequence

- order:

  Multivariate Markov chain order. Default is 2.

- Norm:

  Norm to be used. Default is 2.

## Value

an hommc object

## References

W.-K. Ching et al. / Linear Algebra and its Applications

## Author

Giorgio Spedicato, Deepak Yadav

## Examples

``` r
data <- matrix(c('2', '1', '3', '3', '4', '3', '2', '1', '3', '3', '2', '1', 
               c('2', '4', '4', '4', '4', '2', '3', '3', '1', '4', '3', '3')), 
               ncol = 2, byrow = FALSE)
               
fitHighOrderMultivarMC(data, order = 2, Norm = 2)                
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
#> Lambda2(1,1) : 6.559805e-05
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
#> Lambda2(1,2) : 6.559841e-05
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
#> Lambda2(2,2) : 8.013877e-09
#> P2(2,2) : 
#>   1   2   3   4
#> 1 0 0.0 0.5 0.0
#> 2 0 0.0 0.0 0.2
#> 3 1 0.5 0.0 0.4
#> 4 0 0.5 0.5 0.4
#> 
```
