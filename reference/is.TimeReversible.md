# Check whether a continuous time Markov chain is time reversible

A continuous time Markov chain with stationary distribution \\\pi\\ is
time reversible if it satisfies the detailed balance equations \\\pi_i
Q\_{ij} = \pi_j Q\_{ji}\\ for all states \\i, j\\.

## Usage

``` r
is.TimeReversible(ctmc)
```

## Arguments

- ctmc:

  A [`ctmc`](ctmc-class.md) object.

## Value

`TRUE` if the chain is time reversible, `FALSE` otherwise.

## Details

The equations are checked up to a relative tolerance of
\\\sqrt{\epsilon}\\ (about \\1.5 \times 10^{-8}\\) of the largest flow
\\\pi_i \|Q\_{ij}\|\\, because \\\pi\\ is computed numerically. When the
chain has several stationary distributions they must all satisfy them.

## References

Dobrow, R. P. (2016). Introduction to Stochastic Processes with R.
Wiley.

## See also

Other continuous time Markov chains:
[`ExpectedTime()`](ExpectedTime.md), [`ctmcFit()`](ctmcFit.md),
[`freq2Generator()`](freq2Generator.md),
[`generatorToTransitionMatrix()`](generatorToTransitionMatrix.md),
[`impreciseProbabilityatT()`](impreciseProbabilityatT.md),
[`is.CTMCirreducible()`](is.CTMCirreducible.md),
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
is.TimeReversible(molecularCTMC)
#> [1] TRUE

# a birth-death chain is reversible, a cycle is not
s <- c("a", "b", "c")
birthDeath <- matrix(c(-0.3, 0.3, 0,
                       0.1, -0.5, 0.4,
                       0, 0.7, -0.7), nrow = 3, byrow = TRUE,
                     dimnames = list(s, s))
is.TimeReversible(new("ctmc", states = s, byrow = TRUE, generator = birthDeath))
#> [1] TRUE
cycle <- matrix(c(-1, 1, 0, 0, -1, 1, 1, 0, -1), nrow = 3, byrow = TRUE,
                dimnames = list(s, s))
is.TimeReversible(new("ctmc", states = s, byrow = TRUE, generator = cycle))
#> [1] FALSE
```
