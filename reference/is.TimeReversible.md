# checks if ctmc object is time reversible

The function returns checks if provided function is time reversible

## Usage

``` r
is.TimeReversible(ctmc)
```

## Arguments

- ctmc:

  a ctmc-class object

## Value

Returns a boolean value stating whether ctmc object is time reversible

a boolean value as described above

## References

INTRODUCTION TO STOCHASTIC PROCESSES WITH R, ROBERT P. DOBROW, Wiley

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
is.TimeReversible(molecularCTMC)
#> [1] TRUE
```
