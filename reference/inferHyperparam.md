# Hyperparameters for Bayesian estimation of a Markov chain

The Bayesian estimation of the package uses independent Dirichlet priors
for the rows of the transition matrix. This function infers their
hyperparameters either from a prior transition matrix, taken as the
prior mean, and a scale for each row (its total hyperparameter), or from
a data set, as the transition counts plus one.

## Usage

``` r
inferHyperparam(transMatr = matrix(), scale = numeric(), data = character())
```

## Arguments

- transMatr:

  A prior transition matrix, with the states as dimnames.

- scale:

  A numeric vector with the scale (sum of the hyperparameters) of each
  row of `transMatr`, in the same order.

- data:

  A sequence of states (or a data set accepted by
  [`markovchainFit()`](markovchainFit.md)) from which the
  hyperparameters are inferred.

## Value

A list with element `scaledInference` (from `transMatr` and `scale`) or
`dataInference` (from `data`), the matrix of hyperparameters. Its rows
and columns are sorted alphanumerically by state name.

## Details

Either `data`, or `transMatr` and `scale`, must be given.

## References

Yalamanchi, S. B. and Spedicato, G. A. (2015). Bayesian inference of
first order Markov chains. Vignette of the markovchain package.

## See also

[`markovchainFit()`](markovchainFit.md) (`method = "map"`), which uses
them.

Other Bayesian inference:
[`predictiveDistribution()`](predictiveDistribution.md),
[`priorDistribution()`](priorDistribution.md)

## Author

Sai Bhargav Yalamanchi, Giorgio Spedicato

## Examples

``` r
data(rain, package = "markovchain")
inferHyperparam(data = rain$rain)
#> $dataInference
#>       0 1-5  6+
#> 0   363 127  61
#> 1-5 137  91  69
#> 6+   51  80 125
#> 

weatherStates <- c("sunny", "cloudy", "rain")
weatherMatrix <- matrix(c(0.7, 0.2, 0.1,
                          0.3, 0.4, 0.3,
                          0.2, 0.4, 0.4),
                        byrow = TRUE, nrow = 3,
                        dimnames = list(weatherStates, weatherStates))
inferHyperparam(transMatr = weatherMatrix, scale = c(10, 10, 10))
#> $scaledInference
#>        cloudy rain sunny
#> cloudy      4    3     3
#> rain        4    4     2
#> sunny       2    1     7
#> 
```
