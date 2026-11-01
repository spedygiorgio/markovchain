# Simulate a continuous time Markov chain

Simulates the jumps of a continuous time Markov chain: in each state
\\i\\ the chain waits an exponential time with rate \\-Q\_{ii}\\, then
jumps to \\j \ne i\\ with probability \\-Q\_{ij} / Q\_{ii}\\ (the
embedded jump chain, see
[`generatorToTransitionMatrix()`](generatorToTransitionMatrix.md)).

## Usage

``` r
rctmc(
  n,
  ctmc,
  initDist = numeric(),
  T = 0,
  include.T0 = TRUE,
  out.type = "list"
)
```

## Arguments

- n:

  Number of jumps to simulate: a positive integer, or `Inf` to simulate
  up to time `T`.

- ctmc:

  A [`ctmc`](ctmc-class.md) object.

- initDist:

  Distribution of the initial state, a probability vector in the order
  of `states(ctmc)`. The default, an empty vector, draws the initial
  state uniformly at random.

- T:

  Time horizon: the jumps after time `T` are not returned. The default
  `0` means no horizon, which requires a finite `n`.

- include.T0:

  Logical: include the initial state, at time 0, in the output?

- out.type:

  `"list"` (default) or `"df"`, the form of the output.

## Value

With `out.type = "list"`, an unnamed list of two vectors of equal
length: the states visited (character) and the times at which they were
entered (numeric). With `out.type = "df"`, a data frame with the same
content in the columns `states` and `time`.

## Details

The simulation stops after `n` jumps, at the first jump after `T` (when
`T > 0`), or when an absorbing state (zero exit rate) is reached,
whichever comes first. The generator is used in the orientation given by
the `byrow` slot of `ctmc`.

## References

Anderson, D. F. (2013). Introduction to Stochastic Processes with
Applications in the Biosciences. University of Wisconsin at Madison.

## See also

Other continuous time Markov chains:
[`ExpectedTime()`](ExpectedTime.md), [`ctmcFit()`](ctmcFit.md),
[`freq2Generator()`](freq2Generator.md),
[`generatorToTransitionMatrix()`](generatorToTransitionMatrix.md),
[`impreciseProbabilityatT()`](impreciseProbabilityatT.md),
[`is.CTMCirreducible()`](is.CTMCirreducible.md),
[`is.TimeReversible()`](is.TimeReversible.md),
[`probabilityatT()`](probabilityatT.md),
[`transition2Generator()`](transition2Generator.md)

## Author

Sai Bhargav Yalamanchi

## Examples

``` r
energyStates <- c("sigma", "sigma_star")
gen <- matrix(c(-3, 3, 1, -1), nrow = 2, byrow = TRUE,
              dimnames = list(energyStates, energyStates))
molecularCTMC <- new("ctmc", states = energyStates, byrow = TRUE,
                     generator = gen, name = "Molecular Transition Model")

# all the jumps up to time 1
rctmc(n = Inf, ctmc = molecularCTMC, T = 1)
#> [[1]]
#> [1] "sigma_star" "sigma"      "sigma_star"
#> 
#> [[2]]
#> [1] 0.0000000 0.2886315 0.3588338
#> 
# five jumps from an initial state drawn from c(0.8, 0.2)
rctmc(n = 5, ctmc = molecularCTMC, initDist = c(0.8, 0.2),
      include.T0 = FALSE, out.type = "df")
#>       states              time
#> 1 sigma_star 0.116296115132973
#> 2      sigma 0.632502219907001
#> 3 sigma_star 0.759035449151487
#> 4      sigma  3.45244700975461
#> 5 sigma_star   3.7278060703559
```
