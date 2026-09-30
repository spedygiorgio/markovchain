# Aggregate a Markov chain's state space by Kullback-Leibler minimization

Reduces the state space of a finite, irreducible, aperiodic Markov chain
to `k` macro-states by the spectral-theoretic method of Deng, Mehta and
Meyn (2011): the macro-chain returned is the one whose "lifted" behavior
(each macro-state visit standing in for its micro-states, weighted by
their share of the stationary distribution) is closest, in
Kullback-Leibler divergence rate, to the original chain.

## Usage

``` r
aggregateStates(
  object,
  k = NULL,
  method = c("adaptive", "spectral-bottom-up", "spectral-top-down")
)

# S4 method for class 'markovchain'
aggregateStates(
  object,
  k = NULL,
  method = c("adaptive", "spectral-bottom-up", "spectral-top-down")
)
```

## Arguments

- object:

  A `markovchain` object representing a finite, irreducible, aperiodic
  discrete-time Markov chain with at least 3 states.

- k:

  The number of macro-states to reduce to, an integer between 2 and the
  number of states minus 1. The default, `NULL`, chooses `k`
  automatically via the eigengap heuristic: the transition matrix's
  eigenvalues are sorted by modulus, and `k` is set to the number of
  eigenvalues before the largest relative drop – a standard,
  parameter-free way to guess how many "slow", well-separated modes the
  chain has. This is a genuine automatic selection, unlike PyDTMC's
  `adaptive`, which only picks which of the two algorithms below to run
  for a `k` the caller must still supply.

- method:

  One of `"adaptive"` (the default), `"spectral-bottom-up"` or
  `"spectral-top-down"`. `"spectral-bottom-up"` grows the partition one
  split at a time from a single macro-state, and is the more reliable
  choice for a large reduction (`k` much smaller than the number of
  states); `"spectral-top-down"` starts from every micro-state on its
  own and repeatedly merges the least costly pair, which suits a small
  reduction (`k` close to the number of states). `"adaptive"` follows
  the same rule of thumb as PyDTMC: top-down below 30 states, otherwise
  bottom-up when `k` is at most 30% of the number of states and top-down
  otherwise.

## Value

A named list:

- `partition`:

  A named list of character vectors giving the original state names
  belonging to each macro-state, suitable for passing to
  [`lump`](lump.md) or [`is.lumpable`](is.lumpable.md). Unlike PyDTMC,
  which only labels the reduced chain's states generically (e.g.
  `"ASBU1"`), this traces every macro-state back to the original states
  it stands for.

- `aggregatedChain`:

  The reduced `markovchain` object, row-stochastic, with states named
  after `partition`.

- `klDivergence`:

  The Kullback-Leibler divergence rate (in bits) between `object` and
  the lifted `aggregatedChain`. It is zero (up to rounding) when every
  source state splits its outgoing probability among a destination
  macro-state's members in the same proportions – those of the
  stationary distribution restricted to that macro-state – regardless of
  the source; this is *stronger* than the Kemeny-Snell strong
  lumpability checked by [`is.lumpable`](is.lumpable.md), which only
  requires the macro-to-macro *totals* to agree across sources (see
  Details).

- `method`:

  The method actually used, after resolving `"adaptive"`.

- `k`:

  The number of macro-states actually used, after resolving an automatic
  `NULL`.

## Details

This targets the same problem as [`autoLump`](autoLump.md), but by a
different and more principled route: `autoLump` clusters the leading
eigenvectors with k-means (a generic, randomized heuristic, fixed here
to a deterministic seed only for reproducibility), while
`aggregateStates` greedily minimizes the actual information-theoretic
quantity that measures how much the aggregation distorts the chain's
dynamics.

**When is the divergence exactly zero?** Not simply whenever `object` is
strongly lumpable with respect to `partition` in the sense of
[`is.lumpable`](is.lumpable.md). Strong lumpability only requires that,
for every pair of macro-states, all micro-states in the same source
macro-state have the same *total* probability of moving to the
destination macro-state; it says nothing about how that total is split
among the destination macro-state's own members. The Kullback-Leibler
divergence used here is sensitive to exactly that split: it is zero (up
to rounding) only when every source state distributes its outgoing
probability across a destination macro-state's members in the same
proportions – those of the stationary distribution restricted to that
macro-state – regardless of which source state it is. This is a
genuinely stronger condition, and a strongly lumpable chain need not
satisfy it: the classical Land of Oz weather chain (see the package
vignette), lumped into `Bad_Weather = {rainy, snowy}` and
`Nice_Weather = {nice}`, is strongly lumpable, and `aggregateStates`
correctly recovers that exact partition as optimal and reproduces
[`lump`](lump.md)'s aggregated transition matrix, but its divergence is
strictly positive, because `rainy` and `snowy` split their probability
between `rainy` and `snowy` themselves differently from one another.

Both methods require the chain to be irreducible (for a unique, strictly
positive stationary distribution) and aperiodic (the internal averaging
step used by `"spectral-top-down"` to re-estimate a working stationary
distribution as macro-states are merged assumes convergence, which is
not guaranteed for a periodic chain). Use [`lazyChain`](lazyChain.md)
first to remove periodicity if needed.

## References

Deng, K., Mehta, P. G. and Meyn, S. P. (2011). Optimal Kullback-Leibler
Aggregation via Spectral Theory of Markov Chains. *IEEE Transactions on
Automatic Control*, 56(12).
[doi:10.1109/TAC.2011.2141350](https://doi.org/10.1109/TAC.2011.2141350)

Fill, J. A. (1991). Eigenvalue bounds on convergence to stationarity for
nonreversible Markov chains, with an application to the exclusion
process. *The Annals of Applied Probability*, 1(1).
[doi:10.1214/aoap/1177005981](https://doi.org/10.1214/aoap/1177005981)

## See also

[`autoLump`](autoLump.md), [`lump`](lump.md),
[`is.lumpable`](is.lumpable.md),
[`closestReversible`](closestReversible.md)

## Examples

``` r
# A chain aggregated into two macro-states {a,b}/{c,d} where, in addition
# to being strongly lumpable, "a" and "b" also split their probability
# *within* each destination block identically (0.1/0.1 and 0.4/0.4): this
# stronger property is what makes the divergence exactly zero (see
# Details for a lumpable-but-nonzero counterexample).
statesNames <- c("a", "b", "c", "d")
P <- matrix(c(0.1, 0.1, 0.4, 0.4,
              0.1, 0.1, 0.4, 0.4,
              0.3, 0.3, 0.2, 0.2,
              0.3, 0.3, 0.2, 0.2), byrow = TRUE, nrow = 4,
            dimnames = list(statesNames, statesNames))
mc <- new("markovchain", states = statesNames, transitionMatrix = P)

result <- aggregateStates(mc, k = 2)
result$partition
#> $Macro_1
#> [1] "a" "b"
#> 
#> $Macro_2
#> [1] "c" "d"
#> 
result$klDivergence # zero up to rounding
#> [1] -1.110223e-16
```
