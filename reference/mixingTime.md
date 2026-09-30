# Mixing time of a Markov chain

Estimates the total-variation mixing time of a finite, irreducible,
aperiodic (i.e. ergodic) discrete-time Markov chain: the smallest number
of steps after which the chain's distribution is within `epsilon` of its
stationary distribution, from every possible starting state.

## Usage

``` r
mixingTime(object, epsilon = 0.25, maxIter = 10000L)

# S4 method for class 'markovchain'
mixingTime(object, epsilon = 0.25, maxIter = 10000L)
```

## Arguments

- object:

  A `markovchain` object representing a finite, irreducible, aperiodic
  discrete-time Markov chain.

- epsilon:

  A single number strictly between `0` and `1`: the total-variation
  threshold that counts as "mixed". The classical default `0.25` follows
  Levin and Peres (2017); it is a conventional choice, not a universal
  constant.

- maxIter:

  A single positive integer: the largest \\t\\ that will be tried before
  giving up. This is a safety limit, not a modelling parameter: it
  exists so that a chain which (numerically) mixes only extremely slowly
  reports a clear error instead of looping for an unbounded number of
  iterations.

## Value

A single positive integer, the estimated mixing time
\\t\_{\mathrm{mix}}(\varepsilon)\\. For the trivial one-state chain, `0`
is returned (it is its own stationary distribution).

## Details

For a row-stochastic transition matrix \\P\\ with stationary
distribution \\\pi\\, define the worst-case total variation distance
after \\t\\ steps as \$\$d(t) = \max_i \tfrac{1}{2}\sum_j \|P^t\_{ij} -
\pi_j\|.\$\$ The mixing time returned is
\$\$t\_{\mathrm{mix}}(\varepsilon) = \min\\t \ge 1 : d(t) \le
\varepsilon\\.\$\$

Unlike [`slem`](slem.md), [`spectralGap`](spectralGap.md) and
[`impliedTimescales`](impliedTimescales.md), `mixingTime()` requires
aperiodicity in addition to irreducibility. This is not an arbitrary
restriction carried over from another implementation: for a periodic
chain, \\P^t(i, \cdot)\\ never converges to \\\pi\\ at all (it keeps
cycling through a fixed set of distributions), so \\d(t)\\ does not go
to zero and "the number of steps until \\d(t) \le \varepsilon\\" is
simply undefined for small enough \\\varepsilon\\. [`slem()`](slem.md)
and [`spectralGap()`](spectralGap.md) remain meaningful for periodic
chains because they summarise the transition matrix's spectrum directly,
without reference to a limit that may not exist.

The implementation repeatedly forms \\P^{t+1} = P^t P\\ and checks
\\d(t)\\ after each multiplication, starting from \\t=1\\, until the
threshold is met or `maxIter` is reached. Its time complexity is
\\O(t\_{\mathrm{mix}} \cdot n^3)\\ and its memory use is \\O(n^2)\\ for
a dense \\n\\-state transition matrix: this is a direct, easy-to-audit
computation, not an asymptotically optimal one (a repeated-squaring
scheme would reach a single large power of \\P\\ faster, but would not
let every intermediate \\t\\ be checked against `epsilon` along the
way). It supports both row- and column-stochastic storage.

## References

Levin, D. A. and Peres, Y. (2017). *Markov Chains and Mixing Times*, 2nd
edition. American Mathematical Society.

## See also

[`slem`](slem.md), [`spectralGap`](spectralGap.md),
[`impliedTimescales`](impliedTimescales.md),
[`period`](structuralAnalysis.md), [`is.irreducible`](is.irreducible.md)

## Examples

``` r
statesNames <- c("a", "b")
mc <- new("markovchain",
  states = statesNames,
  transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
    byrow = TRUE, nrow = 2,
    dimnames = list(statesNames, statesNames)))
mixingTime(mc)
#> [1] 3
mixingTime(mc, epsilon = 0.01) # a tighter threshold needs more steps
#> [1] 9
```
