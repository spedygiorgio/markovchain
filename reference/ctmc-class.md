# Continuous time Markov chain class

The S4 class of continuous time Markov chains, described by their
generator (rate) matrix.

## Slots

- `states`:

  Names of the states, equal to the row and column names of the
  generator.

- `byrow`:

  Logical: are the rates out of each state in its row (rows sum to zero)
  or in its column?

- `generator`:

  The square generator matrix: non-negative off-diagonal rates, with
  rows (or columns) summing to zero.

- `name`:

  Optional name of the chain.

## Note

The validity method checks that the rows (or columns) of the generator
sum to zero, rounded to the fifth decimal. If the generator has no state
names, the states are named "s1", "s2", ...

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
plot(molecularCTMC)

```
