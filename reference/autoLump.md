# Automatically aggregate a Markov chain by spectral clustering

Finds an approximate partition by clustering the leading right
eigenvectors of the transition matrix, then returns the forced lumping
over that partition. This is a heuristic for approximate
lumping/metastable aggregation, not a proof of exact lumpability.

## Usage

``` r
autoLump(object, k)

# S4 method for class 'markovchain'
autoLump(object, k)
```

## Arguments

- object:

  A `markovchain` object.

- k:

  Number of macro-states to discover.

## Value

A list with `partition` and `lumped_chain`.
