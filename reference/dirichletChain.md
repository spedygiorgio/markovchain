# Markov chain from a Dirichlet process

Generates a Markov chain whose rows are drawn from a truncated Dirichlet
process with the stick-breaking (GEM) construction, as
`MarkovChain.dirichlet_process()` of PyDTMC does.

## Usage

``` r
dirichletChain(
  n,
  diffusion,
  states = NULL,
  diagonalBias = NULL,
  shiftConcentration = FALSE,
  byrow = TRUE,
  seed = NULL,
  name = "Dirichlet process chain"
)
```

## Arguments

- n:

  The number of states, at least 2. It can be omitted when `states` is
  given.

- diffusion:

  The concentration parameter \\\alpha \> 0\\ of the Dirichlet process.
  Small values concentrate the probability of each row on its first
  states; large values spread it more evenly.

- states:

  An optional character vector of `n` state names. Defaults to
  `as.character(1:n)`.

- diagonalBias:

  An optional positive number \\\beta\\. When given, a draw from
  \\\mathrm{Beta}(\beta, 1)\\ is added to each diagonal entry before the
  row is renormalised, which makes the chain more likely to stay where
  it is; larger values give a stronger bias.

- shiftConcentration:

  If `TRUE`, the columns are reversed, so that the probability
  concentrates on the last states instead of the first ones.

- byrow:

  Whether the transition matrix of the result is stored by rows (the
  default) or by columns.

- seed:

  An optional whole number, as in
  [`randomMarkovChain`](randomMarkovChain.md).

- name:

  The `name` slot of the result.

## Value

A `markovchain` object with `n` states.

## Details

For each row, \\b_1, \ldots, b_n\\ are independent \\\mathrm{Beta}(1,
\alpha)\\ draws and the weights are \$\$w_j = b_j \prod\_{k \< j} (1 -
b_k),\$\$ normalised to sum to one (the truncation at \\n\\ states
leaves out the mass \\\prod_k (1 - b_k)\\). PyDTMC only accepts whole
values of \\\alpha\\ between 1 and \\n\\; any positive value is accepted
here, since the construction is defined for every \\\alpha \> 0\\.

## References

Sethuraman, J. (1994). A constructive definition of Dirichlet priors.
*Statistica Sinica*, 4(2), 639-650.

## See also

[`randomMarkovChain`](randomMarkovChain.md)

## Examples

``` r
dirichletChain(5, diffusion = 2, seed = 1)
#> Dirichlet process chain 
#>  A  5 - dimensional discrete Markov Chain defined by the following states: 
#>  1, 2, 3, 4, 5 
#>  The transition matrix  (by rows)  is defined as follows: 
#>           1         2           3          4          5
#> 1 0.6297278 0.1236413 0.009434792 0.07335242 0.16384368
#> 2 0.2641065 0.1507858 0.453149400 0.01874959 0.11320873
#> 3 0.6825133 0.0346931 0.161828110 0.10215682 0.01880864
#> 4 0.1936219 0.2395930 0.123421772 0.13778836 0.30557497
#> 5 0.1578606 0.3053509 0.218171241 0.29528612 0.02333114
#> 
# a chain that tends to stay in its current state
dirichletChain(5, diffusion = 2, diagonalBias = 5, seed = 1)
#> Dirichlet process chain 
#>  A  5 - dimensional discrete Markov Chain defined by the following states: 
#>  1, 2, 3, 4, 5 
#>  The transition matrix  (by rows)  is defined as follows: 
#>            1          2           3           4          5
#> 1 0.80587070 0.06482367 0.004946549 0.038457799 0.08590129
#> 2 0.13330163 0.57137879 0.228716666 0.009463421 0.05713949
#> 3 0.37722940 0.01917509 0.536737127 0.056462719 0.01039566
#> 4 0.10173955 0.12589530 0.064852562 0.546946757 0.16056583
#> 5 0.08668598 0.16767726 0.119804333 0.162150415 0.46368201
#> 
```
