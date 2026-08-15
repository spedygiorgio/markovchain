# markovchainListFit

Given a data frame or a matrix (rows are observations, by cols the
temporal sequence), it fits a non - homogeneous discrete time markov
chain process (storing row). In particular a markovchainList of size =
ncol - 1 is obtained estimating transitions from the n samples given by
consecutive column pairs.

## Usage

``` r
markovchainListFit(data, byrow = TRUE, laplacian = 0, name)
```

## Arguments

- data:

  Either a matrix or a data.frame or a list object.

- byrow:

  Indicates whether distinc stochastic processes trajectiories are shown
  in distinct rows.

- laplacian:

  Laplacian correction (default 0).

- name:

  Optional name.

## Value

A list containing two slots: estimate (the estimate) name

## Details

If `data` contains `NAs` then the transitions containing `NA` will be
ignored.

## Examples

``` r
# using holson dataset
data(holson)
# fitting a single markovchain
singleMc <- markovchainFit(data = holson[,2:12])
# fitting a markovchainList
mclistFit <- markovchainListFit(data = holson[, 2:12], name = "holsonMcList")
```
