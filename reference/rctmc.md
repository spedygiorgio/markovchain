# rctmc

The function generates random CTMC transitions as per the provided
generator matrix.

## Usage

``` r
rctmc(n, ctmc, initDist = numeric(), T = 0, include.T0 = TRUE,
  out.type = "list")
```

## Arguments

- n:

  The number of samples to generate.

- ctmc:

  The CTMC S4 object.

- initDist:

  The initial distribution of states.

- T:

  The time up to which the simulation runs (all transitions after time T
  are not returned).

- include.T0:

  Flag to determine if start state is to be included.

- out.type:

  "list" or "df"

## Value

Based on out.type, a list or a data frame is returned. The returned list
has two elements - a character vector (states) and a numeric vector
(indicating time of transitions). The data frame is similarly
structured.

## Details

In order to use the T0 argument, set n to Inf.

## References

Introduction to Stochastic Processes with Applications in the
Biosciences (2013), David F. Anderson, University of Wisconsin at
Madison

## See also

[`generatorToTransitionMatrix`](generatorToTransitionMatrix.md),[`ctmc-class`](ctmc-class.md)

## Author

Sai Bhargav Yalamanchi

## Examples

``` r
energyStates <- c("sigma", "sigma_star")
byRow <- TRUE
gen <- matrix(data = c(-3, 3, 1, -1), nrow = 2,
             byrow = byRow, dimnames = list(energyStates, energyStates))
molecularCTMC <- new("ctmc", states = energyStates, 
                     byrow = byRow, generator = gen, 
                     name = "Molecular Transition Model")   

statesDist <- c(0.8, 0.2)
rctmc(n = Inf, ctmc = molecularCTMC, T = 1)
#> [[1]]
#> [1] "sigma_star"
#> 
#> [[2]]
#> [1] 0
#> 
rctmc(n = 5, ctmc = molecularCTMC, initDist = statesDist, include.T0 = FALSE)
#> [[1]]
#> [1] "sigma"      "sigma_star" "sigma"      "sigma_star" "sigma"     
#> 
#> [[2]]
#> [1] 1.413688 1.942383 2.351227 2.805303 5.178264
#> 
```
