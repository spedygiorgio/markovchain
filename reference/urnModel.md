# Build an Ehrenfest urn model Markov chain

Constructs the Ehrenfest diffusion model: `balls` balls are split
between two urns, A and B; at each step, one of the `balls` balls is
chosen uniformly at random and moved to the other urn. The chain tracks
the number of balls in urn A.

## Usage

``` r
urnModel(balls, states = NULL)
```

## Arguments

- balls:

  A single positive integer, the total number of balls. The chain has
  `balls + 1` states, \\0,1,\ldots,\code{balls}\\ (the possible counts
  of balls in urn A).

- states:

  An optional character vector of `balls + 1` state names, in increasing
  order of ball count. Defaults to `as.character(0:balls)`.

## Value

A new, row-stochastic `markovchain` object with `balls + 1` states. From
state \\i\\ (\\0\<i\<\code{balls}\\), \$\$P\_{i,i-1} = i/\code{balls},
\qquad P\_{i,i+1} = 1 - i/\code{balls},\$\$ the probability that the
ball moved was one of the \\i\\ currently in urn A (decreasing A's
count) versus one of the \\\code{balls}-i\\ currently in urn B
(increasing it). States \\0\\ and `balls` (all balls in one urn) are
*reflecting*: the next ball moved must come from the only non-empty urn,
so \\P\_{0,1}=P\_{\code{balls}, \code{balls}-1}=1\\ exactly.

## Details

The Ehrenfest model is the classical example of a chain whose
equilibrium behaviour matches thermodynamic intuition despite every
individual transition being fully reversible: its stationary
distribution is \\\mathrm{Binomial}(\code{balls}, 1/2)\\ (each ball is,
at equilibrium, independently in urn A or B with probability \\1/2\\),
sharply concentrated around \\\code{balls}/2\\ for large `balls` even
though the chain only ever moves one ball at a time and is reflecting,
not absorbing, at the boundaries. It is irreducible and reversible for
every `balls`, but periodic with period \\2\\ (the parity of the ball
count in urn A alternates every step): pass the result through
[`lazyChain`](lazyChain.md) first if an aperiodic chain is needed, e.g.
for [`mixingTime`](mixingTime.md).

## References

Ehrenfest, P. and Ehrenfest, T. (1907). Uber zwei bekannte Einwande
gegen das Boltzmannsche H-Theorem. *Physikalische Zeitschrift*, 8,
311-314.

## See also

[`birthDeath`](birthDeath.md), [`lazyChain`](lazyChain.md)

## Examples

``` r
ehrenfest <- urnModel(balls = 4)
ehrenfest
#> Ehrenfest Urn Model (balls = 4) 
#>  A  5 - dimensional discrete Markov Chain defined by the following states: 
#>  0, 1, 2, 3, 4 
#>  The transition matrix  (by rows)  is defined as follows: 
#>      0   1    2   3    4
#> 0 0.00 1.0 0.00 0.0 0.00
#> 1 0.25 0.0 0.75 0.0 0.00
#> 2 0.00 0.5 0.00 0.5 0.00
#> 3 0.00 0.0 0.75 0.0 0.25
#> 4 0.00 0.0 0.00 1.0 0.00
#> 
steadyStates(ehrenfest) # approximately Binomial(4, 0.5): 1/16 6/16 ...
#>           0    1     2    3      4
#> [1,] 0.0625 0.25 0.375 0.25 0.0625
dbinom(0:4, 4, 0.5)
#> [1] 0.0625 0.2500 0.3750 0.2500 0.0625
```
