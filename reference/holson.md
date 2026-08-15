# Holson data set

A data set containing 1000 life histories trajectories and a categorical
status (1,2,3) observed on eleven evenly spaced steps.

## Usage

``` r
data(holson)
```

## Format

A data frame with 1000 observations on the following 12 variables.

- `id`:

  unique id

- `time1`:

  observed status at i-th time

- `time2`:

  observed status at i-th time

- `time3`:

  observed status at i-th time

- `time4`:

  observed status at i-th time

- `time5`:

  observed status at i-th time

- `time6`:

  observed status at i-th time

- `time7`:

  observed status at i-th time

- `time8`:

  observed status at i-th time

- `time9`:

  observed status at i-th time

- `time10`:

  observed status at i-th time

- `time11`:

  observed status at i-th time

## Source

Private communications

## Details

The example can be used to fit a `markovchain` or a `markovchainList`
object.

## References

Private communications

## Examples

``` r
data(holson)
head(holson)
#>   id time1 time2 time3 time4 time5 time6 time7 time8 time9 time10 time11
#> 1  1     1     1     1     1     1     1     1     1     1      1      1
#> 2  2     1     1     1     1     1     1     1     1     1      1      1
#> 3  3     1     1     1     1     1     1     1     1     1      1      1
#> 4  4     1     1     1     1     2     1     1     1     1      1      1
#> 5  5     1     1     1     1     1     1     1     1     1      1      1
#> 6  6     1     1     1     1     1     1     1     1     1      1      1
```
