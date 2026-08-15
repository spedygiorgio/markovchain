# Continuous time Markov Chains class

The S4 class that describes `ctmc` (continuous time Markov chain)
objects.

## Arguments

- states:

  Name of the states. Must be the same of `colnames` and `rownames` of
  the generator matrix

- byrow:

  TRUE or FALSE. Indicates whether the given matrix is stochastic by
  rows or by columns

- generator:

  Square generator matrix

- name:

  Optional character name of the Markov chain

## Note

1.  `ctmc` classes are written using S4 classes

2.  Validation method is used to assess whether either columns or rows
    totals to zero. Rounding is used up to 5th decimal. If state names
    are not properly defined for a generator `matrix`, coercing to
    `ctmc` object leads to overriding states name with artificial "s1",
    "s2", ... sequence

## Methods

- dim:

  `signature(x = "ctmc")`: method to get the size

- initialize:

  `signature(.Object = "ctmc")`: initialize method

- states:

  `signature(object = "ctmc")`: states method.

- steadyStates:

  `signature(object = "ctmc")`: method to get the steady state vector.

- plot:

  `signature(x = "ctmc", y = "missing")`: plot method for `ctmc` objects

## References

Introduction to Stochastic Processes with Applications in the
Biosciences (2013), David F. Anderson, University of Wisconsin at
Madison. Sai Bhargav Yalamanchi, Giorgio Spedicato

## See also

[`generatorToTransitionMatrix`](generatorToTransitionMatrix.md),[`rctmc`](rctmc.md)

## Examples

``` r
energyStates <- c("sigma", "sigma_star")
byRow <- TRUE
gen <- matrix(data = c(-3, 3,
                       1, -1), nrow = 2,
              byrow = byRow, dimnames = list(energyStates, energyStates))
molecularCTMC <- new("ctmc", states = energyStates, 
                     byrow = byRow, generator = gen, 
                     name = "Molecular Transition Model")
                     steadyStates(molecularCTMC)
#>      sigma sigma_star
#> [1,]  0.25       0.75
if (FALSE) plot(molecularCTMC) # \dontrun{}
```
