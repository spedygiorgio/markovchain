# Returns a generator matrix corresponding to frequency matrix

The function provides interface to calculate generator matrix
corresponding to a frequency matrix and time taken

## Usage

``` r
freq2Generator(P, t = 1, method = "QO", logmethod = "Eigen")
```

## Arguments

- P:

  relative frequency matrix

- t:

  (default value = 1)

- method:

  one among "QO"(Quasi optimaisation), "WA"(weighted adjustment),
  "DA"(diagonal adjustment)

- logmethod:

  method for computation of matrx algorithm (by default : Eigen)

## Value

returns a generator matix with same dimnames

## References

E. Kreinin and M. Sidelnikova: Regularization Algorithms for Transition
Matrices. Algo Research Quarterly 4(1):23-40, 2001

## Examples

``` r
sample <- matrix(c(150,2,1,1,1,200,2,1,2,1,175,1,1,1,1,150),nrow = 4,byrow = TRUE)
sample_rel = rbind((sample/rowSums(sample))[1:dim(sample)[1]-1,],c(rep(0,dim(sample)[1]-1),1)) 
freq2Generator(sample_rel,1)
#>              [,1]        [,2]         [,3] [,4]
#> [1,] -0.024212164  0.01544797  0.008764198    0
#> [2,]  0.006594821 -0.01822834  0.011633520    0
#> [3,]  0.013302567  0.00749703 -0.020799597    0
#> [4,]  0.000000000  0.00000000  0.000000000    0

data(tm_abs)
tm_rel=rbind((tm_abs/rowSums(tm_abs))[1:7,],c(rep(0,7),1))
## Derive quasi optimization generator matrix estimate
freq2Generator(tm_rel,1)
#>              AAA           AA            A          BBB           BB
#> AAA -0.109688198  0.104742772  0.004945426  0.000000000  0.000000000
#> AA   0.006375885 -0.095416840  0.088027218  0.001013738  0.000000000
#> A    0.000000000  0.037605369 -0.139127796  0.092863518  0.002082790
#> BBB  0.000000000  0.003101629  0.043766820 -0.100963213  0.044471253
#> BB   0.000000000  0.004024583  0.000000000  0.043976931 -0.142486437
#> B    0.000000000  0.005844572  0.003289642  0.005803757  0.058923108
#> C    0.000000000  0.000000000  0.000000000  0.000000000  0.006651241
#>      0.000000000  0.000000000  0.000000000  0.000000000  0.000000000
#>                 B            C           D
#> AAA  0.000000e+00  0.000000000 0.000000000
#> AA   0.000000e+00  0.000000000 0.000000000
#> A    1.064427e-05  0.004562575 0.002002900
#> BBB  4.257673e-03  0.001871779 0.003494059
#> BB   8.610403e-02  0.008380896 0.000000000
#> B   -1.932222e-01  0.064440303 0.054920850
#> C    1.547477e-01 -0.362361432 0.200962499
#>      0.000000e+00  0.000000000 0.000000000
```
