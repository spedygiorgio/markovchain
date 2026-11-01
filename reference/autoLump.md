# Aggregate a Markov chain by spectral clustering

Finds a partition of the states into `k` blocks by k-means clustering of
the states on the right eigenvectors of the `k` eigenvalues of largest
modulus, then aggregates the chain over it with
[`lump(force = TRUE)`](lump.md). States that the slow modes of the chain
do not distinguish (metastable groups) end up in the same block.

## Usage

``` r
autoLump(object, k)

# S4 method for class 'markovchain'
autoLump(object, k)
```

## Arguments

- object:

  A [`markovchain`](markovchain-class.md) object.

- k:

  Number of blocks, an integer between 2 and the number of states minus
  1.

## Value

A list with `partition`, a named list of the states of each block
(`"Macro_1"`, `"Macro_2"`, ...), and `lumped_chain`, the aggregated
[`markovchain`](markovchain-class.md) object.

## Details

This is a heuristic for approximate (metastable) aggregation, not a test
of exact lumpability: check the result with
[`is.lumpable()`](is.lumpable.md), or use
[`aggregateStates()`](aggregateStates.md) for the aggregation that
minimizes the Kullback-Leibler divergence rate. The k-means clustering
is run with a fixed seed and 10 starts, so the result is reproducible;
the random number stream of the session is left unchanged.

## See also

Other state aggregation: [`is.lumpable()`](is.lumpable.md),
[`lump()`](lump.md)

## Examples

``` r
# two groups of states that rarely communicate
statesNames <- c("a", "b", "c", "d")
P <- matrix(c(0.6, 0.38, 0.01, 0.01,
              0.3, 0.68, 0.01, 0.01,
              0.01, 0.01, 0.5, 0.48,
              0.01, 0.01, 0.4, 0.58), nrow = 4, byrow = TRUE,
            dimnames = list(statesNames, statesNames))
mc <- new("markovchain", states = statesNames, transitionMatrix = P)
autoLump(mc, k = 2)
#> $partition
#> $partition$Macro_1
#> [1] "a" "b"
#> 
#> $partition$Macro_2
#> [1] "c" "d"
#> 
#> 
#> $lumped_chain
#> Unnamed Markov chain (Lumped) 
#>  A  2 - dimensional discrete Markov Chain defined by the following states: 
#>  Macro_1, Macro_2 
#>  The transition matrix  (by rows)  is defined as follows: 
#>         Macro_1 Macro_2
#> Macro_1    0.98    0.02
#> Macro_2    0.02    0.98
#> 
#> 
```
