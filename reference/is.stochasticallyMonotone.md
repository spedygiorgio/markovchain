# Check whether a Markov chain is stochastically monotone

A Markov chain on ordered states is stochastically monotone if starting
from a higher state never makes the next state stochastically smaller:
for every pair of consecutive states \\i\\, \\i + 1\\ and every \\j\\,
\\\sum\_{l \le j} P\_{i+1, l} \le \sum\_{l \le j} P\_{i, l}\\.

## Usage

``` r
is.stochasticallyMonotone(object)

# S4 method for class 'markovchain'
is.stochasticallyMonotone(object)

# S4 method for class 'matrix'
is.stochasticallyMonotone(object)

# S4 method for class 'ANY'
is.stochasticallyMonotone(object)
```

## Arguments

- object:

  A [`markovchain`](markovchain-class.md) object, or a row-stochastic
  transition matrix.

## Value

`TRUE` if the chain is stochastically monotone, `FALSE` otherwise.

## Details

The states are taken in the order of the transition matrix, and the
inequalities are checked up to a tolerance of \\10^{-10}\\.

## References

Keilson, J. and Kester, A. (1977). Monotone matrices and monotone Markov
processes. Stochastic Processes and their Applications 5(3), 231-241.

## Examples

``` r
# a birth-death chain is monotone
P <- matrix(c(0.7, 0.3, 0,
              0.2, 0.5, 0.3,
              0, 0.4, 0.6), nrow = 3, byrow = TRUE,
            dimnames = list(c("low", "mid", "high"), c("low", "mid", "high")))
is.stochasticallyMonotone(new("markovchain", transitionMatrix = P))
#> [1] TRUE
# reversing the order of the rows breaks monotonicity
is.stochasticallyMonotone(P[3:1, ])
#> [1] FALSE
```
