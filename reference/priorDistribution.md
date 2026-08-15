# priorDistribution

Function to evaluate the prior probability of a transition matrix. It is
based on conjugate priors and therefore a Dirichlet distribution is used
to model the transitions of each state.

## Usage

``` r
priorDistribution(transMatr, hyperparam = matrix())
```

## Arguments

- transMatr:

  The transition matrix whose probability is the parameter of interest.

- hyperparam:

  The hyperparam matrix (optional). If not provided, a default value of
  1 is assumed for each and therefore the resulting probability
  distribution is uniform.

## Value

The log of the probabilities for each state is returned in a numeric
vector. Each number in the vector represents the probability (log) of
having a probability transition vector as specified in corresponding the
row of the transition matrix.

## Details

The states (dimnames) of the transition matrix and the hyperparam may be
in any order.

## Note

This function can be used in conjunction with inferHyperparam. For
example, if the user has a prior data set and a prior transition matrix,
he can infer the hyperparameters using inferHyperparam and then compute
the probability of their prior matrix using the inferred hyperparameters
with priorDistribution.

## References

Yalamanchi SB, Spedicato GA (2015). Bayesian Inference of First Order
Markov Chains. R package version 0.2.5

## See also

[`predictiveDistribution`](predictiveDistribution.md),
[`inferHyperparam`](inferHyperparam.md)

## Author

Sai Bhargav Yalamanchi, Giorgio Spedicato

## Examples

``` r
priorDistribution(matrix(c(0.5, 0.5, 0.5, 0.5), 
                  nrow = 2, 
                  dimnames = list(c("a", "b"), c("a", "b"))), 
                  matrix(c(2, 2, 2, 2), 
                  nrow = 2, 
                  dimnames = list(c("a", "b"), c("a", "b"))))
#>         a         b 
#> 0.4054651 0.4054651 
```
