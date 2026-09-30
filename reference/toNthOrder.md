# Return the n-step transition chain

Returns the `markovchain` object whose transition matrix is
\\P^{\code{order}}\\: from any state, its one-step transition
probabilities are the original chain's `order`-step transition
probabilities.

## Usage

``` r
toNthOrder(object, order)
```

## Arguments

- object:

  A `markovchain` object.

- order:

  A single integer of at least `2`.

## Value

A new `markovchain` object on the same states as `object`, with
transition matrix \\P^{\code{order}}\\.

## Details

This is a thin, discoverability-only wrapper around `object ^ order`
(see [`^,markovchain,numeric-method`](markovchain-class.md)), provided
under this name because PyDTMC's equivalent method is called
`to_nth_order()`. It exists so that the operation is easy to find by
that name; it introduces no new computation; the underlying `^` method
is already \\O(n^3\log(\code{order}))\\ via repeated squaring
(`expm::`[`%^%`](https://rdrr.io/pkg/expm/man/matpow.html)), not a naive
`order`-fold product, so there is nothing to improve on algorithmically
here.

## See also

[`toBoundedChain`](toBoundedChain.md), [`lazyChain`](lazyChain.md)

## Examples

``` r
mc <- new("markovchain", states = c("a", "b"),
          transitionMatrix = matrix(c(0.9, 0.1, 0.3, 0.7), byrow = TRUE, nrow = 2))
identical(unclass(toNthOrder(mc, 5)@transitionMatrix), unclass((mc ^ 5)@transitionMatrix))
#> [1] TRUE
```
