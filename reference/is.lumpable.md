# Check exact lumpability of a Markov chain

Checks whether a Markov chain is (strongly) lumpable with respect to a
partition of its states: for every pair of blocks \\A, B\\ of the
partition, all the states of \\A\\ have the same probability of moving
into \\B\\ (Kemeny and Snell 1960, Theorem 6.3.2). The aggregated
process is then a Markov chain, for every initial distribution.

## Usage

``` r
is.lumpable(object, partition, tol = 1e-10)

# S4 method for class 'markovchain'
is.lumpable(object, partition, tol = 1e-10)
```

## Arguments

- object:

  A [`markovchain`](markovchain-class.md) object.

- partition:

  A named list of character vectors, the blocks of the partition
  (macro-states); every state must belong to exactly one block.

- tol:

  Non-negative tolerance for the equality of the probabilities (default
  `1e-10`).

## Value

`TRUE` if the chain is lumpable with respect to `partition`, `FALSE`
otherwise.

## References

Kemeny, J. G. and Snell, J. L. (1960). Finite Markov Chains. Van
Nostrand.

## See also

Other state aggregation: [`autoLump()`](autoLump.md),
[`lump()`](lump.md)

## Examples

``` r
statesNames <- c("a", "b", "c", "d")
P <- matrix(c(0.1, 0.2, 0.3, 0.4,
              0.2, 0.1, 0.4, 0.3,
              0.3, 0.3, 0.2, 0.2,
              0.2, 0.4, 0.1, 0.3), nrow = 4, byrow = TRUE,
            dimnames = list(statesNames, statesNames))
mc <- new("markovchain", states = statesNames, transitionMatrix = P)
partition <- list(ab = c("a", "b"), cd = c("c", "d"))
is.lumpable(mc, partition) # rows of a and b give 0.3/0.7 to ab/cd
#> [1] TRUE
is.lumpable(mc, list(ac = c("a", "c"), bd = c("b", "d")))
#> [1] FALSE
```
