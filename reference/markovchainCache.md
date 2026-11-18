# Session cache of eigenvalues

Functions of the package that need the eigenvalues of a transition
matrix ([`slem`](slem.md), [`spectralGap`](spectralGap.md),
[`relaxationTime`](relaxationTime.md),
[`impliedTimescales`](impliedTimescales.md),
[`topologicalEntropy`](topologicalEntropy.md),
[`aggregateStates`](aggregateStates.md) with `k = NULL`,
`autoplot(type = "eigenvalues")` and `summary(details = TRUE)`) compute
them once per matrix and keep them for the rest of the R session, so
that calling several of these functions on the same chain does not
repeat the \\O(n^3)\\ eigen-decomposition.

## Details

A stored result is reused only for a matrix bitwise identical to the one
it was computed from (state names are ignored), so results are exactly
the same with or without the cache. The most recently used matrices are
kept within a memory budget.

## Options

- `markovchain.cache`:

  Logical; `FALSE` switches the cache off. Default `TRUE`.

- `markovchain.cacheSize`:

  Memory budget of the cache in bytes, default `64 * 1024^2` (64 MiB). A
  matrix larger than the budget is not stored.

## Examples

``` r
mc <- new("markovchain", states = c("a", "b", "c"),
  transitionMatrix = matrix(c(0.5, 0.3, 0.2,
                              0.1, 0.8, 0.1,
                              0.3, 0.3, 0.4), nrow = 3, byrow = TRUE))
slem(mc)          # computes the eigenvalues
#> [1] 0.5
spectralGap(mc)   # reuses them
#> [1] 0.5

# switch the cache off for the rest of the session
old <- options(markovchain.cache = FALSE)
spectralGap(mc)
#> [1] 0.5
options(old)
```
