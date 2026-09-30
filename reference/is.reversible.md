# Check whether a Markov chain is reversible

Checks whether a finite, irreducible discrete-time Markov chain is
reversible with respect to its (unique) stationary distribution, i.e.
whether it satisfies the detailed balance equations.

## Usage

``` r
is.reversible(object, tolerance = sqrt(.Machine$double.eps))

# S4 method for class 'markovchain'
is.reversible(object, tolerance = sqrt(.Machine$double.eps))
```

## Arguments

- object:

  A `markovchain` object representing a finite, irreducible
  discrete-time Markov chain.

- tolerance:

  A single finite non-negative number. Detailed balance is accepted as
  holding when every pair \\(i,j)\\ satisfies \\\|\pi_i P\_{ij} - \pi_j
  P\_{ji}\| \le \code{tolerance}\\. The default is a small
  numerical-noise tolerance, not a modelling tolerance: it exists to
  absorb floating-point rounding in the eigendecomposition-based
  [`steadyStates`](steadyStates.md) computation, not to declare "almost
  reversible" chains reversible.

## Value

A single logical value, `TRUE` or `FALSE`.

## Details

A chain with transition matrix \\P\\ and stationary distribution \\\pi\\
is reversible if \$\$\pi_i P\_{ij} = \pi_j P\_{ji} \quad \text{for every
} i,j.\$\$ Intuitively, if you started the chain from \\\pi\\ and
watched a long run of it, running the recorded sequence of states
backwards would look statistically identical to running it forwards: at
stationarity, the "flow" of probability from \\i\\ to \\j\\ exactly
balances the flow from \\j\\ to \\i\\.

Only irreducibility is required, not aperiodicity: detailed balance is a
purely algebraic condition on \\P\\ and \\\pi\\ and is perfectly well
defined for periodic chains too. For example, a simple random walk on
any undirected graph (moving to a uniformly random neighbour) is always
reversible, whether or not it happens to be periodic.

A 2-state irreducible chain is always reversible: with only two states,
the single detailed balance equation \\\pi_1 P\_{12} = \pi_2 P\_{21}\\
is just a restatement of the stationarity equation \\\pi P = \pi\\, so
it holds automatically.

Every reversible chain has a real spectrum (all eigenvalues of \\P\\ are
real), which is why [`slem`](slem.md) and
[`spectralGap`](spectralGap.md) are especially easy to interpret for
reversible chains: there are no complex-conjugate eigenvalue pairs to
reason about.

The implementation calls [`steadyStates`](steadyStates.md) once, then
compares the two triangles of the flow matrix \\\pi_i P\_{ij}\\. Its
time complexity is dominated by [`steadyStates()`](steadyStates.md),
plus an additional \\O(n^2)\\ comparison for a dense \\n\\-state
transition matrix. It supports both row- and column-stochastic storage.

## References

Norris, J. R. (1998). *Markov Chains*. Cambridge University Press.

Levin, D. A. and Peres, Y. (2017). *Markov Chains and Mixing Times*, 2nd
edition. American Mathematical Society.

## See also

[`steadyStates`](steadyStates.md),
[`is.irreducible`](is.irreducible.md), [`slem`](slem.md),
[`mixingTime`](mixingTime.md)

## Examples

``` r
# A random walk on a triangle is reversible.
statesNames <- c("a", "b", "c")
triangle <- new("markovchain", states = statesNames,
  transitionMatrix = matrix(c(0, 0.5, 0.5,
                              0.5, 0, 0.5,
                              0.5, 0.5, 0), byrow = TRUE, nrow = 3,
                            dimnames = list(statesNames, statesNames)))
is.reversible(triangle)
#> [1] TRUE

# A directed cycle (states only move "forward") is not reversible: there
# is a net clockwise flow of probability at stationarity.
cycle3 <- new("markovchain", states = statesNames,
  transitionMatrix = matrix(c(0, 1, 0,
                              0, 0, 1,
                              1, 0, 0), byrow = TRUE, nrow = 3,
                            dimnames = list(statesNames, statesNames)))
is.reversible(cycle3)
#> [1] FALSE
```
