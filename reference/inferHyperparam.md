# Function to infer the hyperparameters for Bayesian inference from an a priori matrix or a data set

Since the Bayesian inference approach implemented in the package is
based on conjugate priors, hyperparameters must be provided to model the
prior probability distribution of the chain parameters. The
hyperparameters are inferred from a given a priori matrix under the
assumption that the matrix provided corresponds to the mean (expected)
values of the chain parameters. A scaling factor vector must be provided
too. Alternatively, the hyperparameters can be inferred from a data set.

## Usage

``` r
inferHyperparam(transMatr = matrix(), scale = numeric(), data = character())
```

## Arguments

- transMatr:

  A valid transition matrix, with dimension names.

- scale:

  A vector of scaling factors, each element corresponds to the row names
  of the provided transition matrix transMatr, in the same order.

- data:

  A data set from which the hyperparameters are inferred.

## Value

Returns the hyperparameter matrix in a list.

## Details

transMatr and scale need not be provided if data is provided.

## Note

The hyperparameter matrix returned is such that the row and column names
are sorted alphanumerically, and the elements in the matrix are
correspondingly permuted.

## References

Yalamanchi SB, Spedicato GA (2015). Bayesian Inference of First Order
Markov Chains. R package version 0.2.5

## See also

[`markovchainFit`](markovchainFit.md),
[`predictiveDistribution`](predictiveDistribution.md)

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
weatherMatrix <- matrix(data = c(0.7, 0.2, 0.1, 
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
