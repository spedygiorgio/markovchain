# First passage time distribution

Distribution of the first passage times of a Markov chain started in
`state`: for every state \\j\\ and step \\k\\, the probability that the
chain enters \\j\\ for the first time at step \\k\\.

## Usage

``` r
firstPassage(object, state, n)
```

## Arguments

- object:

  A [`markovchain`](markovchain-class.md) object.

- state:

  Name of the initial state.

- n:

  Number of steps, a positive integer.

## Value

A matrix with `n` rows (steps \\1, \dots, n\\) and one column per state:
entry \\(k, j)\\ is the probability that the first passage to \\j\\
happens at step \\k\\. The column sums are at most one; the remainder is
the probability of a first passage after step `n`, or never.

## Details

For \\j\\ equal to `state` the probabilities are those of the first
return. The algorithm follows the Matlab listings of Feres.

## References

Feres, R. Notes for Math 450: Matlab listings for Markov chains.
Washington University in St. Louis.

## See also

Other first passage and hitting:
[`absorptionProbabilities()`](absorptionProbabilities.md),
[`committorAB()`](committorAB.md),
[`firstPassageMultiple()`](firstPassageMultiple.md),
[`hittingProbabilities()`](hittingProbabilities.md),
[`meanAbsorptionTime()`](meanAbsorptionTime.md),
[`meanFirstPassageTime()`](meanFirstPassageTime.md),
[`meanNumVisits()`](meanNumVisits.md),
[`meanRecurrenceTime()`](meanRecurrenceTime.md),
[`noofVisitsDist()`](noofVisitsDist.md)

## Author

Giorgio Spedicato

## Examples

``` r
simpleMc <- new("markovchain", states = c("a", "b"),
                transitionMatrix = matrix(c(0.4, 0.6, 0.3, 0.7),
                                          nrow = 2, byrow = TRUE))
firstPassage(simpleMc, "b", 20)
#>               a            b
#> 1  0.3000000000 7.000000e-01
#> 2  0.2100000000 1.800000e-01
#> 3  0.1470000000 7.200000e-02
#> 4  0.1029000000 2.880000e-02
#> 5  0.0720300000 1.152000e-02
#> 6  0.0504210000 4.608000e-03
#> 7  0.0352947000 1.843200e-03
#> 8  0.0247062900 7.372800e-04
#> 9  0.0172944030 2.949120e-04
#> 10 0.0121060821 1.179648e-04
#> 11 0.0084742575 4.718592e-05
#> 12 0.0059319802 1.887437e-05
#> 13 0.0041523862 7.549747e-06
#> 14 0.0029066703 3.019899e-06
#> 15 0.0020346692 1.207960e-06
#> 16 0.0014242685 4.831838e-07
#> 17 0.0009969879 1.932735e-07
#> 18 0.0006978915 7.730941e-08
#> 19 0.0004885241 3.092376e-08
#> 20 0.0003419669 1.236951e-08
```
