# Preprogluccacon DNA protein bases sequences

Sequence of bases for preproglucacon DNA protein

## Usage

``` r
data(preproglucacon)
```

## Format

A data frame with 1572 observations on the following 2 variables.

- `V1`:

  a numeric vector, showing original coding

- `preproglucacon`:

  a character vector, showing initial of DNA bases (Adenine, Cytosine,
  Guanine, Thymine)

## Source

Avery Henderson

## References

Averuy Henderson, Fitting markov chain models on discrete time series
such as DNA sequences

## Examples

``` r
data(preproglucacon)
preproglucaconMc<-markovchainFit(data=preproglucacon$preproglucacon)
```
