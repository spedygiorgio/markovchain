# Next-state probabilities and simulation for higher order Markov chains

`higherOrderPredict` returns the distribution of the next state given
the most recent states, and `higherOrderSimulate` draws a sequence of
states, under a higher order model fitted by
[`fitHigherOrder`](fitHigherOrder.md) or [`fitMTD`](fitMTD.md).

## Usage

``` r
higherOrderPredict(fit, history)

higherOrderSimulate(n, fit, t0, include.t0 = FALSE)
```

## Arguments

- fit:

  The list returned by [`fitHigherOrder`](fitHigherOrder.md) or
  [`fitMTD`](fitMTD.md).

- history:

  The most recent states, oldest first: a vector of at least `order`
  states, or a matrix (or data frame) with one history per row.

- n:

  Number of states to simulate.

- t0:

  The states preceding the simulated ones, oldest first; at least
  `order` states.

- include.t0:

  Should `t0` be included at the beginning of the returned sequence?

## Value

`higherOrderPredict`: a named vector of next-state probabilities for a
single history, or a matrix with one row per history.
`higherOrderSimulate`: a character vector of `n` states, preceded by
`t0` if `include.t0 = TRUE`.

## Details

Both functions use the transition probabilities of the fitted model of
order \\k\\, \$\$P(x_t = j \mid x\_{t-1}, \dots, x\_{t-k}) =
\sum\_{i=1}^{k} \lambda_i\\ Q_i\[j, x\_{t-i}\],\$\$ the same that enter
[`higherOrderLogLik`](higherOrderLogLik.md). For a fit returned by
[`fitMTD`](fitMTD.md) all the \\Q_i\\ are the single MTD transition
matrix. Only the last \\k\\ states of a history are used, the last
element being the most recent.

A model of order \\k\\ needs \\k\\ previous states, so `t0` (and every
history) must contain at least `order` states. The probabilities are
normalized to sum to one, which only matters when the weights of a least
squares fit sum to one up to the optimizer's tolerance.

## See also

[`fitHigherOrder`](fitHigherOrder.md), [`fitMTD`](fitMTD.md),
[`higherOrderLogLik`](higherOrderLogLik.md),
[`rmarkovchain`](rmarkovchain.md)

## Examples

``` r
wind <- read.csv(system.file("extdata", "koeberg_wind.csv",
                             package = "markovchain"))$state
fit <- fitMTD(wind, order = 2)
# next wind direction after directions 1 and then 2
higherOrderPredict(fit, c(1, 2))
#>          1          2          3          4 
#> 0.22684612 0.70110514 0.04968651 0.02236223 
# several histories at once
higherOrderPredict(fit, rbind(c(1, 1), c(2, 2), c(4, 1)))
#>               1          2           3            4
#> [1,] 0.83058641 0.06876926 0.007657052 9.298728e-02
#> [2,] 0.03568195 0.90132362 0.062994430 2.053392e-56
#> [3,] 0.64957465 0.05223115 0.018490373 2.797038e-01
# simulate one day of hourly directions starting from the last two observed
set.seed(1)
higherOrderSimulate(24, fit, t0 = tail(wind, 2))
#>  [1] "1" "1" "1" "4" "4" "1" "2" "2" "2" "2" "2" "2" "2" "2" "2" "2" "2" "1" "1"
#> [20] "1" "2" "2" "2" "2"

# the same works for fitHigherOrder()
data(rain)
fit2 <- fitHigherOrder(rain$rain, order = 2, method = "mle")
higherOrderPredict(fit2, c("0", "6+"))
#>         0       1-5        6+ 
#> 0.2464501 0.3040798 0.4494701 
```
