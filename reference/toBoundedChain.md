# Apply a boundary condition to a Markov chain's first and last state

Replaces the transition rows of a `markovchain` object's first and last
state (in `states(object)` order) with an absorbing, reflecting, or
semi-reflecting rule, leaving every other row unchanged.

## Usage

``` r
toBoundedChain(object, boundaryCondition)

# S4 method for class 'markovchain'
toBoundedChain(object, boundaryCondition)
```

## Arguments

- object:

  A `markovchain` object with at least 2 states.

- boundaryCondition:

  Either:

  - the string `"absorbing"`: the first and last state each become
    absorbing (\\P\_{11}=1\\, \\P\_{nn}=1\\);

  - the string `"reflecting"`: the first state moves to the second with
    certainty and the last state moves to the second-to-last with
    certainty (\\P\_{12}=1\\, \\P\_{n,n-1}=1\\);

  - a single number \\\beta\in\[0,1\]\\, the *semi-reflecting* case: the
    first state stays with probability \\1-\beta\\ and moves to the
    second state with probability \\\beta\\ (\\P\_{11}=1-\beta\\,
    \\P\_{12}=\beta\\), and symmetrically the last state stays with
    probability \\1-\beta\\ and moves to the second-to-last with
    probability \\\beta\\. \\\beta=0\\ is the absorbing case and
    \\\beta=1\\ is the reflecting case.

## Value

A new `markovchain` object, row-stochastic, on the same states as
`object`, identical to `object` except in its first and last transition
rows.

## Details

This function assumes – as is standard for a boundary condition – that
`states(object)` is meaningfully ordered along a line, first state to
last state, as it would be e.g. for [`birthDeath`](birthDeath.md) or any
other chain built to represent a bounded random walk. It does not check
this (there is no general way to check it from the transition matrix
alone) and applies the same first/last-row replacement regardless of
`object`'s actual structure; only the two boundary rows are ever
touched, so applying it to a chain whose states are not linearly ordered
simply reinterprets whichever states happen to be listed first and last.

Unlike [`gamblersRuin`](gamblersRuin.md), which is absorbing at both
ends by construction and cannot be un-done, `toBoundedChain()` can be
applied to any existing chain and with any of the three conditions,
including reflecting or semi-reflecting ones that
[`gamblersRuin()`](gamblersRuin.md) does not offer directly.

The implementation touches only 2 of the \\n\\ rows and is \\O(n)\\ time
and memory beyond copying the transition matrix.

## See also

[`birthDeath`](birthDeath.md), [`gamblersRuin`](gamblersRuin.md)

## Examples

``` r
bd <- birthDeath(p = c(0.3, 0.4, 0.5), q = c(0.2, 0.3, 0.1))

absorbed <- toBoundedChain(bd, "absorbing")
absorbed@transitionMatrix[1, ]
#> 1 2 3 4 
#> 1 0 0 0 
absorbed@transitionMatrix[4, ]
#> 1 2 3 4 
#> 0 0 0 1 

reflected <- toBoundedChain(bd, "reflecting")
reflected@transitionMatrix[1, ]
#> 1 2 3 4 
#> 0 1 0 0 

semiReflected <- toBoundedChain(bd, 0.25)
semiReflected@transitionMatrix[1, ]
#>    1    2    3    4 
#> 0.75 0.25 0.00 0.00 
```
