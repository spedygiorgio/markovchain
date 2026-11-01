# Fit a non-homogeneous Markov chain

Estimates a non-homogeneous discrete time Markov chain, with one
transition matrix for each pair of consecutive times, from a panel of
trajectories observed at the same times.

## Usage

``` r
markovchainListFit(data, byrow = TRUE, laplacian = 0, name)
```

## Arguments

- data:

  The trajectories: a matrix or data frame with one trajectory per row
  and one time per column (or the transpose, see `byrow`), or a list of
  trajectories (character vectors), possibly of different lengths.

- byrow:

  Logical: for a matrix or data frame, is each trajectory a row
  (default) or a column?

- laplacian:

  Laplace smoothing: a constant added to every transition count before
  normalizing (default 0).

- name:

  Optional name of the fitted
  [`markovchainList`](markovchainList-class.md).

## Value

A list with element `estimate`, the fitted
[`markovchainList`](markovchainList-class.md); an empty list if no
transition can be counted.

## Details

The \\k\\-th transition matrix is the maximum likelihood estimate from
the transitions observed between times \\k\\ and \\k + 1\\ (columns
\\k\\ and \\k + 1\\, or positions \\k\\ and \\k + 1\\ of the
trajectories long enough). Transitions with a missing value (`NA`) are
ignored, and pairs of times without any complete transition are skipped.

## See also

[`markovchainFit()`](markovchainFit.md) for a homogeneous chain,
[`rmarkovchain()`](rmarkovchain.md) to simulate from the fit.

## Examples

``` r
data(holson)
# a homogeneous chain for the whole panel
singleMc <- markovchainFit(data = holson[, 2:12])
# one transition matrix for each pair of consecutive years
mclistFit <- markovchainListFit(data = holson[, 2:12], name = "holsonMcList")
mclistFit$estimate
#> holsonMcList  list of Markov chain(s) 
#> Markovchain  1 
#> Unnamed Markov chain 
#>  A  3 - dimensional discrete Markov Chain defined by the following states: 
#>  1, 2, 3 
#>  The transition matrix  (by rows)  is defined as follows: 
#>            1          2         3
#> 1 0.94609164 0.05390836 0.0000000
#> 2 0.26356589 0.62790698 0.1085271
#> 3 0.02325581 0.18604651 0.7906977
#> 
#> Markovchain  2 
#> Unnamed Markov chain 
#>  A  3 - dimensional discrete Markov Chain defined by the following states: 
#>  1, 2, 3 
#>  The transition matrix  (by rows)  is defined as follows: 
#>           1         2         3
#> 1 0.9323410 0.0676590 0.0000000
#> 2 0.2551724 0.5103448 0.2344828
#> 3 0.0000000 0.0862069 0.9137931
#> 
#> Markovchain  3 
#> Unnamed Markov chain 
#>  A  3 - dimensional discrete Markov Chain defined by the following states: 
#>  1, 2, 3 
#>  The transition matrix  (by rows)  is defined as follows: 
#>            1          2           3
#> 1 0.94765840 0.04820937 0.004132231
#> 2 0.26119403 0.66417910 0.074626866
#> 3 0.01428571 0.13571429 0.850000000
#> 
#> Markovchain  4 
#> Unnamed Markov chain 
#>  A  3 - dimensional discrete Markov Chain defined by the following states: 
#>  1, 2, 3 
#>  The transition matrix  (by rows)  is defined as follows: 
#>           1          2           3
#> 1 0.9172414 0.07724138 0.005517241
#> 2 0.1678322 0.60839161 0.223776224
#> 3 0.0000000 0.03030303 0.969696970
#> 
#> Markovchain  5 
#> Unnamed Markov chain 
#>  A  3 - dimensional discrete Markov Chain defined by the following states: 
#>  1, 2, 3 
#>  The transition matrix  (by rows)  is defined as follows: 
#>           1          2         3
#> 1 0.9579100 0.04208999 0.0000000
#> 2 0.1972789 0.70748299 0.0952381
#> 3 0.0000000 0.17073171 0.8292683
#> 
#> Markovchain  6 
#> Unnamed Markov chain 
#>  A  3 - dimensional discrete Markov Chain defined by the following states: 
#>  1, 2, 3 
#>  The transition matrix  (by rows)  is defined as follows: 
#>           1          2           3
#> 1 0.9462990 0.05224964 0.001451379
#> 2 0.1677019 0.61490683 0.217391304
#> 3 0.0000000 0.06666667 0.933333333
#> 
#> Markovchain  7 
#> Unnamed Markov chain 
#>  A  3 - dimensional discrete Markov Chain defined by the following states: 
#>  1, 2, 3 
#>  The transition matrix  (by rows)  is defined as follows: 
#>             1          2           3
#> 1 0.949926362 0.04860088 0.001472754
#> 2 0.193103448 0.71034483 0.096551724
#> 3 0.005681818 0.15909091 0.835227273
#> 
#> Markovchain  8 
#> Unnamed Markov chain 
#>  A  3 - dimensional discrete Markov Chain defined by the following states: 
#>  1, 2, 3 
#>  The transition matrix  (by rows)  is defined as follows: 
#>            1          2         3
#> 1 0.92581602 0.07418398 0.0000000
#> 2 0.06707317 0.73780488 0.1951220
#> 3 0.00000000 0.04320988 0.9567901
#> 
#> Markovchain  9 
#> Unnamed Markov chain 
#>  A  3 - dimensional discrete Markov Chain defined by the following states: 
#>  1, 2, 3 
#>  The transition matrix  (by rows)  is defined as follows: 
#>           1          2          3
#> 1 0.9637795 0.03622047 0.00000000
#> 2 0.2247191 0.70786517 0.06741573
#> 3 0.0000000 0.17647059 0.82352941
#> 
#> Markovchain  10 
#> Unnamed Markov chain 
#>  A  3 - dimensional discrete Markov Chain defined by the following states: 
#>  1, 2, 3 
#>  The transition matrix  (by rows)  is defined as follows: 
#>           1          2         3
#> 1 0.9585890 0.04141104 0.0000000
#> 2 0.1318681 0.74725275 0.1208791
#> 3 0.0000000 0.06626506 0.9337349
#> 
```
