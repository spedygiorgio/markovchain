# Aggregate a Markov chain over a partition

Coarsens a Markov chain to a reduced state space. By default the
function requires exact lumpability. With `force = TRUE`, it performs an
approximate aggregation using stationary weights when available and
arithmetic averages for macro-states with zero stationary mass.

## Usage

``` r
lump(object, partition, force = FALSE)

# S4 method for class 'markovchain'
lump(object, partition, force = FALSE)
```

## Arguments

- object:

  A `markovchain` object.

- partition:

  A named list of character vectors defining macro-states.

- force:

  If `FALSE`, stop unless the chain is exactly lumpable. If `TRUE`,
  return a weighted approximate lumping.

## Value

A `markovchain` object on the macro-state space.
