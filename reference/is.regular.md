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
elements strictly positive

## References

Matrix Analysis. Roger A.Horn, Charles R.Johnson. 2nd edition. Corollary
8.5.8, Theorem 8.5.9

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
