# Bayesian predictive probability of new data

Log-probability of observing the sequence `newData` after the sequence
`stringchar`, under a first-order Markov chain whose rows have
independent Dirichlet priors: the likelihood of `newData` averaged over
the posterior given `stringchar`.

## Usage

``` r
predictiveDistribution(stringchar, newData, hyperparam = matrix())
```

## Arguments

- stringchar:

  The observed sequence of states on which the posterior is based (a
  character vector).

- newData:

  The sequence whose predictive probability is computed.

- hyperparam:

  Matrix of Dirichlet hyperparameters, one row per state, with the
  states as dimnames (in any order). The default, an empty matrix, sets
  all of them to 1 (uniform priors).

## Value

A number, the logarithm of the predictive probability.

## Details

With conjugate Dirichlet priors the predictive probability has a closed
form, a ratio of products of gamma functions (Strelioff et al. 2007),
which is what is computed.

## References

Strelioff, C. C., Crutchfield, J. P. and Hubler, A. W. (2007). Inferring
Markov chains: Bayesian estimation, model comparison, entropy rate, and
out-of-class modeling. Physical Review E 76, 011106.
[doi:10.1103/PhysRevE.76.011106](https://doi.org/10.1103/PhysRevE.76.011106)

Yalamanchi, S. B. and Spedicato, G. A. (2015). Bayesian inference of
first order Markov chains. Vignette of the markovchain package.

## See also

Other Bayesian inference: [`inferHyperparam()`](inferHyperparam.md),
[`priorDistribution()`](priorDistribution.md)

## Author

Sai Bhargav Yalamanchi

## Examples

``` r
sequence <- c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a",
              "a", "b", "b", "b", "a")
hyperMatrix <- matrix(c(1, 2, 1, 4), nrow = 2,
                      dimnames = list(c("a", "b"), c("a", "b")))
predProb <- predictiveDistribution(sequence[1:10], sequence[11:17],
                                   hyperparam = hyperMatrix)
predProb
#> [1] -4.489759
# the order of the states in the hyperparameter matrix does not matter
hyperMatrix2 <- hyperMatrix[c(2, 1), c(2, 1)]
predictiveDistribution(sequence[1:10], sequence[11:17],
                       hyperparam = hyperMatrix2) == predProb
#> [1] TRUE
```
