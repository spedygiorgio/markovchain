# Random Markov chain

Generates a Markov chain whose transition probabilities are drawn at
random, optionally with a given number of zero probabilities and with
some probabilities fixed in advance.

## Usage

``` r
randomMarkovChain(
  n,
  states = NULL,
  zeros = 0L,
  mask = NULL,
  byrow = TRUE,
  seed = NULL,
  name = "Random Markov chain"
)
```

## Arguments

- n:

  The number of states. It can be omitted when `states` is given.

- states:

  An optional character vector of `n` state names. Defaults to
  `as.character(1:n)`.

- zeros:

  The number of transition probabilities, among those not fixed by
  `mask`, that are set to zero. Every row whose probabilities are not
  all fixed keeps at least one positive free entry, which bounds `zeros`
  from above; a larger value is an error.

- mask:

  An optional `n x n` matrix of fixed transition probabilities: `NA`
  marks the entries to draw at random, any other value (in \\\[0, 1\]\\)
  is kept as it is. In each row the fixed values must not sum to more
  than one. A row whose fixed values sum to one gets zero in its `NA`
  entries; in any other row the free entries share the remaining
  probability. With `byrow = FALSE` the mask is read by columns, like
  the transition matrix.

- byrow:

  Whether the transition matrix of the result is stored by rows (the
  default) or by columns.

- seed:

  An optional whole number. When given, the chain is generated after
  `set.seed(seed)` and the caller's random number stream is restored
  afterwards, so the result is reproducible without affecting later
  random draws; when `NULL` (the default), the current stream is used,
  so [`set.seed()`](https://rdrr.io/r/base/Random.html) beforehand works
  as usual.

- name:

  The `name` slot of the result.

## Value

A `markovchain` object with `n` states.

## Details

The algorithm follows `MarkovChain.random()` of PyDTMC. In each row not
completely fixed by `mask`, one free entry, chosen at random, is
reserved to be positive. The `zeros` zero entries are then chosen at
random among the remaining free entries of the whole matrix, the other
free entries are drawn from a uniform distribution on \\(0, 1)\\, and
the free entries of each row are rescaled so that, together with the
fixed ones, they sum to one. The rows are normalised uniforms, which is
not the uniform distribution on the simplex; use
[`dirichletChain`](dirichletChain.md) or draw the rows yourself if that
matters.

## See also

[`dirichletChain`](dirichletChain.md),
[`identityChain`](identityChain.md)

## Examples

``` r
randomMarkovChain(4, seed = 1)
#> Random Markov chain 
#>  A  4 - dimensional discrete Markov Chain defined by the following states: 
#>  1, 2, 3, 4 
#>  The transition matrix  (by rows)  is defined as follows: 
#>            1          2         3         4
#> 1 0.09022035 0.28142773 0.3073326 0.3210193
#> 2 0.38455404 0.02644750 0.1644149 0.4245836
#> 3 0.41063439 0.08953367 0.3346371 0.1651949
#> 4 0.31280384 0.08357720 0.2355974 0.3680216
#> 
# 6 of the 16 transition probabilities are zero
sum(randomMarkovChain(4, zeros = 6, seed = 1)@transitionMatrix == 0)
#> [1] 6
# state "b" moves to "a" with probability 0.5; the rest is random
m <- matrix(NA, 3, 3)
m[2, 1] <- 0.5
randomMarkovChain(states = c("a", "b", "c"), mask = m, seed = 1)
#> Random Markov chain 
#>  A  3 - dimensional discrete Markov Chain defined by the following states: 
#>  a, b, c 
#>  The transition matrix  (by rows)  is defined as follows: 
#>           a         b         c
#> a 0.3728717 0.3688408 0.2582876
#> b 0.5000000 0.4693052 0.0306948
#> c 0.1887605 0.6184614 0.1927781
#> 
```
