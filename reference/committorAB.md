# Committor probabilities

Probability that a Markov chain hits the set of states `A` before the
set `B` (the committor function), from every initial state or from one.

## Usage

``` r
committorAB(object, A, B, p = 1)
```

## Arguments

- object:

  A [`markovchain`](markovchain-class.md) object.

- A, B:

  Disjoint sets of states, given as indices (positions in
  `states(object)`).

- p:

  Optional index of the initial state. If omitted, the probabilities for
  all initial states are returned.

## Value

A numeric vector with the probability for each initial state, or a
single number if `p` is given.

## Details

The committor \\q\\ solves the linear system \\q_i = \sum_j P\_{ij}
q_j\\ for the states outside \\A \cup B\\, with \\q = 1\\ on `A` and \\q
= 0\\ on `B`. The system has a unique solution when `A` or `B` can be
reached from every state.

## See also

Other first passage and hitting:
[`absorptionProbabilities()`](absorptionProbabilities.md),
[`firstPassage()`](firstPassage.md),
[`firstPassageMultiple()`](firstPassageMultiple.md),
[`hittingProbabilities()`](hittingProbabilities.md),
[`meanAbsorptionTime()`](meanAbsorptionTime.md),
[`meanFirstPassageTime()`](meanFirstPassageTime.md),
[`meanNumVisits()`](meanNumVisits.md),
[`meanRecurrenceTime()`](meanRecurrenceTime.md),
[`noofVisitsDist()`](noofVisitsDist.md)

## Examples

``` r
transMatr <- matrix(c(0, 0, 0, 1, 0.5,
                      0.5, 0, 0, 0, 0,
                      0.5, 0, 0, 0, 0,
                      0, 0.2, 0.4, 0, 0,
                      0, 0.8, 0.6, 0, 0.5), nrow = 5)
object <- new("markovchain", states = c("a", "b", "c", "d", "e"),
              transitionMatrix = transMatr)
# probability of reaching "e" before "c", from each state
committorAB(object, 5, 3)
#>         a         b         c         d         e 
#> 0.4444444 0.8888889 0.0000000 0.4444444 1.0000000 
committorAB(object, 5, 3, p = 1)
#>         a 
#> 0.4444444 
```
