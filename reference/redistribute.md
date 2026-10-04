# Evolution of a distribution over time

Propagates an initial distribution of the states of a discrete-time
Markov chain forward in time and returns the whole trajectory of
distributions.

## Usage

``` r
redistribute(object, steps, initial = NULL, lastOnly = FALSE)

# S4 method for class 'markovchain'
redistribute(object, steps, initial = NULL, lastOnly = FALSE)
```

## Arguments

- object:

  A `markovchain` object.

- steps:

  A single non-negative whole number: the number of steps to propagate.
  With `steps = 0` only the initial distribution is returned.

- initial:

  The initial distribution. Either `NULL` (default), for the uniform
  distribution over the states; a single state name, for a point mass on
  that state; or a numeric vector of non-negative probabilities summing
  to one. A named numeric vector is matched to the states by name, an
  unnamed one by position.

- lastOnly:

  Logical. If `TRUE`, only the distribution after `steps` steps is
  returned. Defaults to `FALSE`.

## Value

If `lastOnly = FALSE` (default), a numeric matrix with `steps + 1` rows
and one column per state: row \\t\\ (labelled `"t"`, from `"0"`) holds
the distribution after \\t\\ steps, so the first row is the initial
distribution. Otherwise, a named numeric vector with the distribution
after `steps` steps.

## Details

If \\\mu_0\\ is the initial distribution (a row vector) and \\P\\ the
row-stochastic transition matrix, the distribution after \\t\\ steps is
\$\$\mu_t = \mu\_{t-1} P = \mu_0 P^t.\$\$

The implementation propagates the vector step by step, at a cost of
\\O(n^2)\\ per step for a dense \\n\\-state chain, instead of forming
\\P^t\\; this is what makes the whole trajectory available at no extra
cost. Both row- and column-stochastic storage are supported. Each
distribution is renormalized after every step to prevent round-off from
accumulating over long horizons.

This mirrors PyDTMC's `redistribute()`, with the same defaults (uniform
initial distribution, output including the initial one). For chains that
converge, the rows approach the stationary distribution (see
[`steadyStates`](steadyStates.md)); for periodic chains they do not,
which is the expected behaviour and not an error.

## See also

[`steadyStates`](steadyStates.md), [`mixingTime`](mixingTime.md),
[`autoplot.markovchain`](autoplot.markovchain.md)

## Examples

``` r
statesNames <- c("a", "b")
mc <- new("markovchain",
  states = statesNames,
  transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
    byrow = TRUE, nrow = 2,
    dimnames = list(statesNames, statesNames)))
redistribute(mc, steps = 5, initial = "a")
#>         a       b
#> 0 1.00000 0.00000
#> 1 0.70000 0.30000
#> 2 0.52000 0.48000
#> 3 0.41200 0.58800
#> 4 0.34720 0.65280
#> 5 0.30832 0.69168
redistribute(mc, steps = 50, initial = c(a = 0.2, b = 0.8), lastOnly = TRUE)
#>    a    b 
#> 0.25 0.75 
```
