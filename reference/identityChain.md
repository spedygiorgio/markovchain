# Identity Markov chain

Builds the Markov chain whose transition matrix is the identity: every
state is absorbing, so the chain never moves.

## Usage

``` r
identityChain(n, states = NULL, name = "Identity chain")
```

## Arguments

- n:

  The number of states. It can be omitted when `states` is given.

- states:

  An optional character vector of `n` state names. Defaults to
  `as.character(1:n)`.

- name:

  The `name` slot of the result.

## Value

A `markovchain` object with `n` states.

## Details

The chain is the neutral element of the product of transition matrices
(`identityChain(n) * mc` equals `mc` for a chain `mc` on the same
states) and the extreme case of [`lazyChain`](lazyChain.md). Every state
is its own closed class, so every distribution is stationary.

## See also

[`randomMarkovChain`](randomMarkovChain.md), [`lazyChain`](lazyChain.md)

## Examples

``` r
identityChain(3)
#> Identity chain 
#>  A  3 - dimensional discrete Markov Chain defined by the following states: 
#>  1, 2, 3 
#>  The transition matrix  (by rows)  is defined as follows: 
#>   1 2 3
#> 1 1 0 0
#> 2 0 1 0
#> 3 0 0 1
#> 
identityChain(states = c("a", "b"))
#> Identity chain 
#>  A  2 - dimensional discrete Markov Chain defined by the following states: 
#>  a, b 
#>  The transition matrix  (by rows)  is defined as follows: 
#>   a b
#> a 1 0
#> b 0 1
#> 
absorbingStates(identityChain(3))
#> [1] "1" "2" "3"
```
