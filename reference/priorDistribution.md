# Prior density of a transition matrix

Log-density of each row of a transition matrix under independent
Dirichlet priors, the conjugate priors used by the Bayesian estimation
of the package.

## Usage

``` r
priorDistribution(transMatr, hyperparam = matrix())
```

## Arguments

- transMatr:

  The transition matrix, with the states as dimnames.

- hyperparam:

  Matrix of Dirichlet hyperparameters, with the states as dimnames (in
  any order). The default, an empty matrix, sets all of them to 1, the
  uniform prior.

## Value

A named numeric vector with, for each state, the logarithm of the
Dirichlet density of its row of `transMatr`.

## Note

With a prior data set and a prior transition matrix, the hyperparameters
can be inferred with [`inferHyperparam()`](inferHyperparam.md) and the
prior density of the matrix then evaluated with these hyperparameters.

## References

Yalamanchi, S. B. and Spedicato, G. A. (2015). Bayesian inference of
first order Markov chains. Vignette of the markovchain package.

## See also

Other Bayesian inference: [`inferHyperparam()`](inferHyperparam.md),
[`predictiveDistribution()`](predictiveDistribution.md)

## Author

Sai Bhargav Yalamanchi, Giorgio Spedicato

## Examples

``` r
priorDistribution(matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2,
                         dimnames = list(c("a", "b"), c("a", "b"))),
                  matrix(c(2, 2, 2, 2), nrow = 2,
                         dimnames = list(c("a", "b"), c("a", "b"))))
#>         a         b 
#> 0.4054651 0.4054651 
```
