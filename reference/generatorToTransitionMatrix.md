# Transition matrix of the embedded jump chain

Computes, from the generator \\Q\\ of a continuous time Markov chain,
the transition matrix of its embedded (jump) chain: the probability of
jumping from \\i\\ to \\j \ne i\\ is \\-Q\_{ij} / Q\_{ii}\\.

## Usage

``` r
generatorToTransitionMatrix(gen, byrow = TRUE)
```

## Arguments

- gen:

  A generator matrix: non-negative off-diagonal rates, with rows (or
  columns, if `byrow = FALSE`) summing to zero.

- byrow:

  Logical: are the rates out of each state in its row (default) or in
  its column? The result has the same orientation.

## Value

The transition matrix of the embedded chain, with the dimnames of `gen`.
A state with zero exit rate (absorbing) gets probability 1 of staying
where it is.

## References

Anderson, D. F. (2013). Introduction to Stochastic Processes with
Applications in the Biosciences. University of Wisconsin at Madison.

## See also

Other continuous time Markov chains:
[`ExpectedTime()`](ExpectedTime.md), [`ctmcFit()`](ctmcFit.md),
[`freq2Generator()`](freq2Generator.md),
[`impreciseProbabilityatT()`](impreciseProbabilityatT.md),
[`is.CTMCirreducible()`](is.CTMCirreducible.md),
[`is.TimeReversible()`](is.TimeReversible.md),
[`probabilityatT()`](probabilityatT.md), [`rctmc()`](rctmc.md),
[`transition2Generator()`](transition2Generator.md)

## Author

Sai Bhargav Yalamanchi

## Examples

``` r
energyStates <- c("sigma", "sigma_star")
gen <- matrix(c(-3, 3, 1, -1), nrow = 2, byrow = TRUE,
              dimnames = list(energyStates, energyStates))
generatorToTransitionMatrix(gen)
#>            sigma sigma_star
#> sigma          0          1
#> sigma_star     1          0
```
