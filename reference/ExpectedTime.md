# Returns expected hitting time from state i to state j

Returns expected hitting time from state i to state j

## Usage

``` r
ExpectedTime(C,i,j,useRCpp)
```

## Arguments

- C:

  A CTMC S4 object

- i:

  Initial state i

- j:

  Final state j

- useRCpp:

  logical whether to use Rcpp

## Value

A numerical value that returns expected hitting times from i to j

## Details

According to the theorem, holding times for all states except j should
be greater than 0.

## References

Markovchains, J. R. Norris, Cambridge University Press

## Author

Vandit Jain

## Examples

``` r
states <- c("a","b","c","d")
byRow <- TRUE
gen <- matrix(data = c(-1, 1/2, 1/2, 0, 1/4, -1/2, 0, 1/4, 1/6, 0, -1/3, 1/6, 0, 0, 0, 0),
nrow = 4,byrow = byRow, dimnames = list(states,states))
ctmc <- new("ctmc",states = states, byrow = byRow, generator = gen, name = "testctmc")
ExpectedTime(ctmc,1,4,TRUE)
#> [1] 7
```
