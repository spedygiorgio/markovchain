# Check whether a continuous time Markov chain is irreducible

A continuous time Markov chain is irreducible if every state can be
reached from every other state, which is the case if and only if its
embedded jump chain is irreducible.

## Usage

``` r
is.CTMCirreducible(ctmc)
```

## Arguments

- ctmc:

  A [`ctmc`](ctmc-class.md) object.

## Value

`TRUE` if the chain is irreducible, `FALSE` otherwise.

## References

Sigman, K. (2009). Continuous-time Markov chains. Lecture notes,
Columbia University.

## See also

Other continuous time Markov chains:
[`ExpectedTime()`](ExpectedTime.md), [`ctmcFit()`](ctmcFit.md),
[`freq2Generator()`](freq2Generator.md),
[`generatorToTransitionMatrix()`](generatorToTransitionMatrix.md),
[`impreciseProbabilityatT()`](impreciseProbabilityatT.md),
[`is.TimeReversible()`](is.TimeReversible.md),
[`probabilityatT()`](probabilityatT.md), [`rctmc()`](rctmc.md),
[`transition2Generator()`](transition2Generator.md)

## Author

Vandit Jain

## Examples

``` r
energyStates <- c("sigma", "sigma_star")
gen <- matrix(c(-3, 3, 1, -1), nrow = 2, byrow = TRUE,
              dimnames = list(energyStates, energyStates))
molecularCTMC <- new("ctmc", states = energyStates, byrow = TRUE,
                     generator = gen, name = "Molecular Transition Model")
is.CTMCirreducible(molecularCTMC)
#> [1] TRUE
```
