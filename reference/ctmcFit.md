# Fit a continuous time Markov chain

Maximum likelihood estimation of a continuous time Markov chain from one
observed trajectory: the sequence of states visited and the times at
which they were entered.

## Arguments

- data:

  A list of two vectors of equal length: the states visited (character)
  and the times at which they were entered (numeric, increasing, the
  first equal to 0).

- byrow:

  Logical: return the generator, and the transition matrix of the
  embedded jump chain in the confidence intervals, with the rates out of
  each state in its row (default) or in its column?

- name:

  Optional name of the fitted chain.

- confidencelevel:

  Confidence level of the intervals (default 0.95).

## Value

A list with two elements:

- `estimate`:

  the fitted [`ctmc`](ctmc-class.md) object;

- `errors`:

  a list with `dtmcConfidenceInterval`, the confidence level and the
  lower and upper endpoint matrices of the transition probabilities of
  the embedded jump chain (as in
  [`markovchainFit()`](markovchainFit.md)), and
  `lambdaConfidenceInterval`, the lower and upper endpoints of the exit
  rate of each state.

## Details

The exit rate of state \\i\\ is estimated as the number of departures
from \\i\\ divided by the total time spent in \\i\\, and the jump
probabilities by the transition frequencies of the sequence of states;
the generator is their product. The interval for each rate is the normal
approximation \\\hat\lambda_i (1 \pm z / \sqrt{n_i})\\, with \\n_i\\ the
number of departures from \\i\\, truncated at 0.

## References

Yalamanchi, S. B. and Spedicato, G. A. (2015). Continuous time Markov
chains. Vignette of the markovchain package.

## See also

Other continuous time Markov chains:
[`ExpectedTime()`](ExpectedTime.md),
[`freq2Generator()`](freq2Generator.md),
[`generatorToTransitionMatrix()`](generatorToTransitionMatrix.md),
[`impreciseProbabilityatT()`](impreciseProbabilityatT.md),
[`is.CTMCirreducible()`](is.CTMCirreducible.md),
[`is.TimeReversible()`](is.TimeReversible.md),
[`probabilityatT()`](probabilityatT.md), [`rctmc()`](rctmc.md),
[`transition2Generator()`](transition2Generator.md)

## Author

Sai Bhargav Yalamanchi

## Examples

``` r
data <- list(c("a", "b", "c", "a", "b", "a", "c", "b", "c"),
             c(0, 0.8, 2.1, 2.4, 4, 5, 5.9, 8.2, 9))
fit <- ctmcFit(data)
fit$estimate
#> An object of class "ctmc"
#> Slot "states":
#> [1] "a" "b" "c"
#> 
#> Slot "byrow":
#> [1] TRUE
#> 
#> Slot "generator":
#>            a          b          c
#> a -0.9090909  0.6060606  0.3030303
#> b  0.3225806 -0.9677419  0.6451613
#> c  0.3846154  0.3846154 -0.7692308
#> 
#> Slot "name":
#> [1] ""
#> 
fit$errors$lambdaConfidenceInterval
#> $lowerEndpointVector
#> [1] 0 0 0
#> 
#> $upperEndpointVector
#> [1] 1.937805 2.062825 1.835311
#> 
```
