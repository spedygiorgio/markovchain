# Sales Demand Sequences

Sales demand sequences of five products (A, B, C, D, E). Each row
corresponds to a sequence. First row corresponds to Sequence A, Second
row to Sequence B and so on.

## Usage

``` r
data("sales")
```

## Format

An object of class `matrix` (inherits from `array`) with 269 rows and 5
columns.

## Details

The example can be used to fit High order multivariate markov chain.

## Examples

``` r
data("sales")
# fitHighOrderMultivarMC(seqMat = sales, order = 2, Norm = 2)
```
