# CD4 cells counts on HIV Infects between zero and six month

This is the table shown in Craig and Sendi paper showing zero and six
month CD4 cells count in six brakets

## Usage

``` r
data(craigsendi)
```

## Format

A 3 x 3 contingency table of counts, with rows and columns named
`"0-49"`, `"50-74"` and `"75-UP"`.

## Source

Estimation of the transition matrix of a discrete time Markov chain,
Bruce A. Craig and Peter P. Sendi, Health Economics 11, 2002.

## Details

Rows represent counts at the beginning, cols represent counts after six
months.

## References

see source

## Examples

``` r
data(craigsendi)
csMc<-as(craigsendi, "markovchain")
steadyStates(csMc)
#>           0-49      50-74      75-UP
#> [1,] 0.8343668 0.07659214 0.08904103
```
