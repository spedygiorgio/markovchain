# Check exact lumpability of a Markov chain

Verifies the strong lumpability condition with respect to a partition of
the state space. For every pair of macro-states, all micro-states in the
same source macro-state must have the same total probability of moving
to the destination macro-state.

## Usage

``` r
is.lumpable(object, partition, tol = 1e-10)

# S4 method for class 'markovchain'
is.lumpable(object, partition, tol = 1e-10)
```

## Arguments

- object:

  A `markovchain` object.

- partition:

  A named list of character vectors defining macro-states.

- tol:

  Non-negative numerical tolerance for equality checks.

## Value

A logical value.

## References

Kemeny, J. G. and Snell, J. L. (1960). *Finite Markov Chains*.
