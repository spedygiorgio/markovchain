# Expected reward collected before hitting a set of states

Expected total reward collected by a Markov chain started in `state` at
the steps \\1, \dots, n\\ it spends outside the set `A` before entering
it, when state \\i\\ pays `rewards[i]` at each visit.

## Usage

``` r
expectedRewardsBeforeHittingA(markovchain, A, state, rewards, n)
```

## Arguments

- markovchain:

  A [`markovchain`](markovchain-class.md) object.

- A:

  Names of the states of the target set.

- state:

  Name of the initial state, not in `A`.

- rewards:

  Numeric vector with the reward of each state, in the order of
  `states(markovchain)` (the rewards of the states in `A` are not used).

- n:

  Number of steps, a non-negative integer.

## Value

A number, the expected reward.

## Details

With \\Q\\ the transition matrix restricted to the states outside `A`,
the value is \\\sum\_{k=1}^{n} (Q^k r)\_{s}\\, with \\s\\ the initial
state: the reward of the initial state at time 0 is not included.

## See also

Other rewards: [`expectedRewards()`](expectedRewards.md)

## Author

Sai Bhargav Yalamanchi, Vandit Jain

## Examples

``` r
transMatr <- matrix(c(0.5, 0.3, 0.2,
                      0.2, 0.5, 0.3,
                      0, 0, 1), nrow = 3, byrow = TRUE)
mc <- new("markovchain", states = c("a", "b", "c"),
          transitionMatrix = transMatr)
# reward 1 per step spent in "a" or "b" before absorption in "c"
expectedRewardsBeforeHittingA(mc, A = "c", state = "a",
                              rewards = c(1, 1, 0), n = 100)
#> [1] 3.210526
```
