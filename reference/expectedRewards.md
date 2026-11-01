# Expected cumulative reward of a Markov chain

Expected total reward collected by a Markov chain over the times \\0, 1,
\dots, n\\, for every initial state, when state \\i\\ pays `rewards[i]`
at each visit.

## Usage

``` r
expectedRewards(markovchain, n, rewards)
```

## Arguments

- markovchain:

  A [`markovchain`](markovchain-class.md) object.

- n:

  Number of steps, a non-negative integer.

- rewards:

  Numeric vector with the reward of each state, in the order of
  `states(markovchain)`.

## Value

A numeric vector with one value per initial state.

## Details

The values solve the dynamic programming recursion \\v_0 = r\\, \\v_k =
r + P v\_{k-1}\\ (Gallager 2013, Section 3.5), so that \\v_n =
\sum\_{k=0}^{n} P^k r\\.

## References

Gallager, R. G. (2013). Stochastic Processes: Theory for Applications.
Cambridge University Press.

## See also

Other rewards:
[`expectedRewardsBeforeHittingA()`](expectedRewardsBeforeHittingA.md)

## Author

Vandit Jain

## Examples

``` r
transMatr <- matrix(c(0.99, 0.01, 0.01, 0.99), nrow = 2, byrow = TRUE)
simpleMc <- new("markovchain", states = c("a", "b"),
                transitionMatrix = transMatr)
expectedRewards(simpleMc, 1, c(0, 1))
#> [1] 0.01 1.99
```
