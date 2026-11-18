# Check if a DTMC is regular

Function to check wether a DTCM is regular

## Usage

``` r
is.regular(object)
```

## Arguments

- object:

  a markovchain object

## Value

A boolean value

## Details

A Markov chain is regular if some of the powers of its matrix has all
elements strictly positive. A finite chain is regular if and only if it
is irreducible and aperiodic (Horn and Johnson, Theorem 8.5.3 and
Corollary 8.5.8), and this is how it is checked, from the communicating
classes and the period ([`period()`](structuralAnalysis.md)); the work
is \\O(n^2)\\ for a dense \\n\\-state chain.

## References

Matrix Analysis. Roger A.Horn, Charles R.Johnson. 2nd edition. Theorem
8.5.3, Corollary 8.5.8, Theorem 8.5.9

## See also

[`is.irreducible`](is.irreducible.md)

## Author

Ignacio Cordón

## Examples

``` r
P <- matrix(c(0.5,  0.25, 0.25,
              0.5,     0, 0.5,
              0.25, 0.25, 0.5), nrow = 3)
colnames(P) <- rownames(P) <- c("R","N","S")
ciao <- as(P, "markovchain")
is.regular(ciao)
#> [1] TRUE
```
