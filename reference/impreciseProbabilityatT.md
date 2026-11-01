# Lower transition probability of an imprecise continuous time Markov chain

For an imprecise continuous time Markov chain, whose transition rates
are only known to lie in intervals, computes the lower probability of
being in state `i` at time `s`, as a function of the state at time `t`
(Krak et al. 2017).

## Usage

``` r
impreciseProbabilityatT(C, i, t = 0, s, error = 10^-3, useRCpp = TRUE)
```

## Arguments

- C:

  An `ictmc` object: its slot `Q` is a generator and row \\k\\ of the
  rates may be scaled by any factor in `range[k, ]`.

- i:

  Index (position in the states of `C`) of the state whose lower
  probability is computed.

- t:

  Initial time (default 0).

- s:

  Final time, larger than `t`.

- error:

  Maximum approximation error (default 0.001).

- useRCpp:

  Logical: use the C++ implementation (default)? It requires `t` and `s`
  to be integer-valued.

## Value

A numeric vector with one element per state \\x\\: the lower probability
that the chain is in state `i` at time `s` given that it is in \\x\\ at
time `t`. With degenerate intervals (`range` equal to 1) it is column
`i` of \\\exp((s - t) Q)\\, up to `error`.

## Details

The lower transition operator is approximated by \\n\\ steps \\g
\leftarrow g + \delta \underline{Q} g\\, with \\\delta = (s - t) / n\\
and \\n\\ chosen so that the error is below `error`, starting from the
indicator of state `i`; \\\underline{Q}\\ takes, in each row, the
smaller of the values given by the two extremes of the interval of
rates.

## References

Krak, T., De Bock, J. and Siebes, A. (2017). Imprecise continuous-time
Markov chains. International Journal of Approximate Reasoning 88,
452-528.

## See also

Other continuous time Markov chains:
[`ExpectedTime()`](ExpectedTime.md), [`ctmcFit()`](ctmcFit.md),
[`freq2Generator()`](freq2Generator.md),
[`generatorToTransitionMatrix()`](generatorToTransitionMatrix.md),
[`is.CTMCirreducible()`](is.CTMCirreducible.md),
[`is.TimeReversible()`](is.TimeReversible.md),
[`probabilityatT()`](probabilityatT.md), [`rctmc()`](rctmc.md),
[`transition2Generator()`](transition2Generator.md)

## Author

Vandit Jain

## Examples

``` r
states <- c("n", "y")
Q <- matrix(c(-1, 1, 1, -1), nrow = 2, byrow = TRUE,
            dimnames = list(states, states))
range <- matrix(c(1/52, 3/52, 1/2, 2), nrow = 2, byrow = TRUE)
ictmc <- new("ictmc", states = states, Q = Q, range = range,
             name = "testictmc")
impreciseProbabilityatT(ictmc, i = 2, t = 0, s = 1, error = 1e-3)
#> [1] 0.008259774 0.140983489
```
