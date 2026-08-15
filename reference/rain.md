# Alofi island daily rainfall

Rainfall measured in Alofi Island

## Usage

``` r
data(rain)
```

## Format

A data frame with 1096 observations on the following 2 variables.

- `V1`:

  a numeric vector, showing original coding

- `rain`:

  a character vector, showing daily rainfall millilitres brackets

## Source

Avery Henderson

## References

Avery Henderson, Fitting markov chain models on discrete time series
such as DNA sequences

## Examples

``` r
data(rain)
rainMc<-markovchainFit(data=rain$rain)
```
