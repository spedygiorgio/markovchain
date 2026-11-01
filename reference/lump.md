# Aggregate a Markov chain over a partition

Builds the chain on the blocks (macro-states) of a partition of the
states. By default the chain must be exactly lumpable with respect to
the partition (see [`is.lumpable()`](is.lumpable.md)); with
`force = TRUE` an approximate aggregation is returned for any partition.

## Usage

``` r
lump(object, partition, force = FALSE)

# S4 method for class 'markovchain'
lump(object, partition, force = FALSE)
```

## Arguments

- object:

  A [`markovchain`](markovchain-class.md) object.

- partition:

  A named list of character vectors, the blocks of the partition; the
  names become the states of the result.

- force:

  Logical: if `FALSE` (default) stop unless the chain is exactly
  lumpable; if `TRUE` aggregate anyway.

## Value

A row-stochastic [`markovchain`](markovchain-class.md) object on the
blocks.

## Details

The probability of moving from block \\A\\ to block \\B\\ is the average
over the states of \\A\\ of their probability of moving into \\B\\,
weighted by the stationary distribution (the average of the stationary
distributions if there are several; equal weights within a block with
zero stationary mass). For a lumpable chain all the states of a block
give the same value, so the weights do not matter; otherwise the result
is the usual stationary-weighted approximation, which preserves the
stationary distribution of the blocks.

## References

Kemeny, J. G. and Snell, J. L. (1960). Finite Markov Chains. Van
Nostrand.

## See also

Other state aggregation: [`autoLump()`](autoLump.md),
[`is.lumpable()`](is.lumpable.md)

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
lump(mc, partition)
#> Unnamed Markov chain (Lumped) 
#>  A  2 - dimensional discrete Markov Chain defined by the following states: 
#>  ab, cd 
#>  The transition matrix  (by rows)  is defined as follows: 
#>     ab  cd
#> ab 0.3 0.7
#> cd 0.6 0.4
#> 
# an approximate aggregation over a partition that is not lumpable
lump(mc, list(ac = c("a", "c"), bd = c("b", "d")), force = TRUE)
#> Unnamed Markov chain (Lumped) 
#>  A  2 - dimensional discrete Markov Chain defined by the following states: 
#>  ac, bd 
#>  The transition matrix  (by rows)  is defined as follows: 
#>           ac        bd
#> ac 0.4543103 0.5456897
#> bd 0.4395833 0.5604167
#> 
```
