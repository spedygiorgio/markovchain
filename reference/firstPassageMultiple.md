# First passage time distribution of a set of states

Distribution of the first passage time into a set of states: the
probability that a Markov chain started in `state` enters `set` for the
first time at step \\k = 1, \dots, n\\.

## Usage

``` r
firstPassageMultiple(object, state, set, n)
```

## Arguments

- object:

  A [`markovchain`](markovchain-class.md) object.

- state:

  Name of the initial state.

- set:

  Names of the states of the target set.

- n:

  Number of steps, a positive integer.

## Value

A one-column matrix (column `"set"`) with `n` rows: row \\k\\ is the
probability that the first passage into `set` happens at step \\k\\.

## References

Feres, R. Notes for Math 450: Matlab listings for Markov chains.
Washington University in St. Louis.

Gallager, R. G. 6.262 Discrete Stochastic Processes, course notes,
chapter 5. MIT OpenCourseWare.

## See also

Other first passage and hitting:
[`absorptionProbabilities()`](absorptionProbabilities.md),
[`committorAB()`](committorAB.md), [`firstPassage()`](firstPassage.md),
[`hittingProbabilities()`](hittingProbabilities.md),
[`meanAbsorptionTime()`](meanAbsorptionTime.md),
[`meanFirstPassageTime()`](meanFirstPassageTime.md),
[`meanNumVisits()`](meanNumVisits.md),
[`meanRecurrenceTime()`](meanRecurrenceTime.md),
[`noofVisitsDist()`](noofVisitsDist.md)

## Author

Vandit Jain

## Examples

``` r
statesNames <- c("a", "b", "c")
markovB <- new("markovchain", states = statesNames, transitionMatrix =
  matrix(c(0.2, 0.5, 0.3,
           0, 1, 0,
           0.1, 0.8, 0.1), nrow = 3, byrow = TRUE,
         dimnames = list(statesNames, statesNames)))
firstPassageMultiple(markovB, "a", c("b", "c"), 4)
#>      set
#> 1 0.8000
#> 2 0.4000
#> 3 0.1190
#> 4 0.0379
```
