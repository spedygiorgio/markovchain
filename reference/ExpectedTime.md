# Expected hitting time of a state in a continuous time Markov chain

Expected time for a continuous time Markov chain started in state `i` to
reach state `j` for the first time.

## Usage

``` r
ExpectedTime(C, i, j, useRCpp = TRUE)
```

## Arguments

- C:

  A [`ctmc`](ctmc-class.md) object.

- i, j:

  Indices (positions in `states(C)`) of the starting and target states.

- useRCpp:

  Logical: solve the linear system in C++ (default) or with
  [`solve()`](https://rdrr.io/r/base/solve.html)?

## Value

A number, the expected hitting time of `j` from `i` (0 if `i == j`).

## Details

The expected hitting times \\m_k\\ of `j` from the states \\k \ne j\\
solve \\\sum\_{l \ne j} Q\_{kl} m_l = -1\\ (Norris 1997, Section 3.3).
Every state other than `j` must have a positive exit rate, otherwise an
error is given; the times are finite only if `j` is reachable from every
state, otherwise the system is singular.

## References

Norris, J. R. (1997). Markov Chains. Cambridge University Press.

## See also

Other continuous time Markov chains: [`ctmcFit()`](ctmcFit.md),
[`freq2Generator()`](freq2Generator.md),
[`generatorToTransitionMatrix()`](generatorToTransitionMatrix.md),
[`impreciseProbabilityatT()`](impreciseProbabilityatT.md),
[`is.CTMCirreducible()`](is.CTMCirreducible.md),
[`is.TimeReversible()`](is.TimeReversible.md),
[`probabilityatT()`](probabilityatT.md), [`rctmc()`](rctmc.md),
[`transition2Generator()`](transition2Generator.md)

## Author

Vandit Jain

## Examples

``` r
states <- c("a", "b", "c", "d")
gen <- matrix(c(-1, 1/2, 1/2, 0,
                1/4, -1/2, 0, 1/4,
                1/6, 0, -1/3, 1/6,
                0, 0, 0, 0), nrow = 4, byrow = TRUE,
              dimnames = list(states, states))
ctmc <- new("ctmc", states = states, byrow = TRUE, generator = gen,
            name = "testctmc")
ExpectedTime(ctmc, 1, 4) # from "a" to the absorbing state "d"
#> [1] 7
```
