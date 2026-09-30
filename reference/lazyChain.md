# Build a lazy version of a Markov chain

Constructs the "lazy" chain associated with a `markovchain` object: at
every step, stay put with probability `alpha` and otherwise take a step
of the original chain.

## Usage

``` r
lazyChain(object, alpha = 0.5)

# S4 method for class 'markovchain'
lazyChain(object, alpha = 0.5)
```

## Arguments

- object:

  A `markovchain` object.

- alpha:

  A single number in \\\[0,1\]\\, the probability of staying in the
  current state at each step. The default, `0.5`, matches the usual
  textbook "lazy random walk" construction.

## Value

A new `markovchain` object with transition matrix \\L = \alpha I +
(1-\alpha) P\\, the same states, and the same row/column-stochastic
storage convention (`byrow`) as `object`.

## Details

For a transition matrix \\P\\ (in whichever storage convention `object`
already uses) and a laziness parameter \\\alpha\in\[0,1\]\\, the lazy
chain's transition matrix is \$\$L = \alpha I + (1-\alpha) P.\$\$

Laziness is a standard device for forcing aperiodicity without changing
where the chain can go or its stationary distribution:

- \\L\\ has the \*same\* stationary distribution as \\P\\ (if \\\pi
  P=\pi\\ then \\\pi L = \alpha\pi + (1-\alpha)\pi P = \pi\\), and the
  same communicating classes, since \\L\_{ij}\>0 \iff P\_{ij}\>0\\ for
  \\i\ne j\\.

- For \\0\<\alpha\<1\\, \\L\\ is aperiodic even if \\P\\ is periodic,
  because \\L\_{ii}=\alpha\>0\\ for every state \\i\\ rules out any
  period greater than \\1\\. This is why [`mixingTime`](mixingTime.md),
  which requires aperiodicity, is often applied to `lazyChain(object)`
  rather than to a periodic `object` directly (see
  [`mixingTime`](mixingTime.md)'s own documentation for why it rejects
  periodic chains outright rather than lazifying them automatically).

- Every non-trivial eigenvalue of \\L\\ is \\\alpha +
  (1-\alpha)\lambda\\ for the corresponding eigenvalue \\\lambda\\ of
  \\P\\: laziness shrinks the whole non-trivial spectrum towards
  \\\alpha\\, so [`slem`](slem.md) and
  [`impliedTimescales`](impliedTimescales.md) generally get \*worse\*
  (mixing gets slower) as `alpha` increases towards `1`.

\\\alpha=0\\ returns \\P\\ unchanged; \\\alpha=1\\ returns the identity
matrix (a chain that never moves).

## References

Levin, D. A. and Peres, Y. (2017). *Markov Chains and Mixing Times*, 2nd
edition. American Mathematical Society.

## See also

[`subchain`](subchain.md), [`mixingTime`](mixingTime.md),
[`slem`](slem.md)

## Examples

``` r
# A 2-cycle is periodic (period 2); its lazy version is aperiodic.
statesNames <- c("a", "b")
cycle2 <- new("markovchain", states = statesNames,
  transitionMatrix = matrix(c(0, 1, 1, 0), byrow = TRUE, nrow = 2,
                            dimnames = list(statesNames, statesNames)))
period(cycle2)
#> [1] 2
lazyCycle2 <- lazyChain(cycle2, alpha = 0.5)
period(lazyCycle2)
#> [1] 1
steadyStates(cycle2)
#>        a   b
#> [1,] 0.5 0.5
steadyStates(lazyCycle2) # unchanged by laziness
#>        a   b
#> [1,] 0.5 0.5
```
