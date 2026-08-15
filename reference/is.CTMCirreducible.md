# Check if CTMC is irreducible

This function verifies whether a CTMC object is irreducible

## Usage

``` r
is.CTMCirreducible(ctmc)
```

## Arguments

- ctmc:

  a ctmc-class object

## Value

a boolean value as described above.

## References

Continuous-Time Markov Chains, Karl Sigman, Columbia University

## Author

Vandit Jain

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
is.CTMCirreducible(molecularCTMC)
#> [1] TRUE
```
