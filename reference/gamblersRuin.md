# Build a gambler's ruin Markov chain

Constructs the classic gambler's ruin chain: a gambler with a fortune
between `0` and `upperBound` wins each round (and gains one unit) with
probability `prob`, otherwise loses one unit; play stops as soon as the
fortune reaches `0` (ruin) or `upperBound` (the gambler's target).

## Usage

``` r
gamblersRuin(upperBound, prob, states = NULL)
```

## Arguments

- upperBound:

  A single positive integer: the fortune at which the gambler stops
  (having won). The chain has `upperBound + 1` states,
  \\0,1,\ldots,\code{upperBound}\\.

- prob:

  A single number in \\\[0,1\]\\: the probability of winning an
  individual round (moving up by one unit) while the fortune is strictly
  between `0` and `upperBound`.

- states:

  An optional character vector of `upperBound + 1` state names, in
  increasing order of fortune. Defaults to `as.character(0:upperBound)`.

## Value

A new, row-stochastic `markovchain` object with `upperBound + 1` states.
States `"0"` and `as.character(upperBound)` are absorbing; every
interior state \\i\\ has \\P\_{i,i+1}=\code{prob}\\ and
\\P\_{i,i-1}=1-\code{prob}\\.

## Details

This is the special case of [`birthDeath`](birthDeath.md) with constant
birth probability `prob` and constant death probability `1-prob` at
every interior state, together forced to be *absorbing* rather than
merely reflecting at the two ends – which is why it is provided as its
own constructor rather than expressed purely in terms of
[`birthDeath()`](birthDeath.md), which cannot produce absorbing
boundaries by itself (see [`toBoundedChain`](toBoundedChain.md) for
turning any chain's ends absorbing or reflecting after construction).

With `prob != 0.5`, the classical ruin probability of reaching `0`
before `upperBound`, starting from fortune \\i\\, is
\$\$P(\text{ruin}\mid X_0=i) =
\frac{\left(\frac{1-\code{prob}}{\code{prob}}\right)^{i} -
\left(\frac{1-\code{prob}}{\code{prob}}\right)^{\code{upperBound}}}
{1-\left(\frac{1-\code{prob}}{\code{prob}}\right)^{\code{upperBound}}},\$\$
and \\i/\code{upperBound}\\ when `prob = 0.5`; this is a standard
textbook result (see Norris (1998), Section 1.3) and is not itself
computed by this function, but can be read off from
[`absorptionProbabilities`](absorptionProbabilities.md) applied to the
returned chain.

## References

Norris, J. R. (1998). *Markov Chains*. Cambridge University Press.

## See also

[`birthDeath`](birthDeath.md),
[`absorptionProbabilities`](absorptionProbabilities.md),
[`toBoundedChain`](toBoundedChain.md)

## Examples

``` r
ruin <- gamblersRuin(upperBound = 5, prob = 0.4)
ruin
#> Gambler's Ruin (upperBound = 5) 
#>  A  6 - dimensional discrete Markov Chain defined by the following states: 
#>  0, 1, 2, 3, 4, 5 
#>  The transition matrix  (by rows)  is defined as follows: 
#>     0   1   2   3   4   5
#> 0 1.0 0.0 0.0 0.0 0.0 0.0
#> 1 0.6 0.0 0.4 0.0 0.0 0.0
#> 2 0.0 0.6 0.0 0.4 0.0 0.0
#> 3 0.0 0.0 0.6 0.0 0.4 0.0
#> 4 0.0 0.0 0.0 0.6 0.0 0.4
#> 5 0.0 0.0 0.0 0.0 0.0 1.0
#> 
absorbingStates(ruin)
#> [1] "0" "5"
```
