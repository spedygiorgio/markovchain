# Transition probabilities of a continuous time Markov chain at time t

Computes the transition matrix \\P(t) = \exp(t Q)\\ of a continuous time
Markov chain with generator \\Q\\, or one of its rows.

## Usage

``` r
probabilityatT(C, t, x0, useRCpp = TRUE)
```

## Arguments

- C:

  A [`ctmc`](ctmc-class.md) object.

- t:

  Time, a non-negative number.

- x0:

  Optional index (position in `states(C)`) of the initial state.

- useRCpp:

  Logical: compute the matrix exponential in C++ (default) or with
  [`expm::expm()`](https://rdrr.io/pkg/expm/man/expm.html)?

## Value

Without `x0`, the matrix \\P(t)\\, with the states as dimnames, whose
entry \\(i, j)\\ is the probability of being in state `j` at time `t`
when starting from `i`. With `x0`, row `x0` of that matrix: the
distribution of the state at time `t` when starting from `x0`.

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
[`is.TimeReversible()`](is.TimeReversible.md), [`rctmc()`](rctmc.md),
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
probabilityatT(ctmc, 1)
#>            a          b         c          d
#> a 0.41546882 0.24714119 0.2703605 0.06702946
#> b 0.12357060 0.63939068 0.0348290 0.20220972
#> c 0.09012017 0.02321933 0.7411205 0.14553997
#> d 0.00000000 0.00000000 0.0000000 1.00000000
probabilityatT(ctmc, 1, x0 = 1)
#> [1] 0.41546882 0.24714119 0.27036052 0.06702946
```
