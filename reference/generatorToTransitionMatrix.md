# Function to obtain the transition matrix from the generator

The transition matrix of the embedded DTMC is inferred from the CTMC's
generator

## Usage

``` r
generatorToTransitionMatrix(gen, byrow = TRUE)
```

## Arguments

- gen:

  The generator matrix

- byrow:

  Flag to determine if rows (columns) sum to 0

## Value

Returns the transition matrix.

## References

Introduction to Stochastic Processes with Applications in the
Biosciences (2013), David F. Anderson, University of Wisconsin at
Madison

## See also

[`rctmc`](rctmc.md),[`ctmc-class`](ctmc-class.md)

## Author

Sai Bhargav Yalamanchi

## Examples

``` r
energyStates <- c("sigma", "sigma_star")
byRow <- TRUE
gen <- matrix(data = c(-3, 3, 1, -1), nrow = 2,
              byrow = byRow, dimnames = list(energyStates, energyStates))
generatorToTransitionMatrix(gen)
#>            sigma sigma_star
#> sigma          0          1
#> sigma_star     1          0
```
