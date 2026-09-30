# Closest reversible approximation of a Markov chain

Finds the Markov chain closest to a given one among those that are
reversible with respect to a fixed stationary distribution.

## Usage

``` r
closestReversible(
  object,
  stationaryDistribution = NULL,
  tolerance = sqrt(.Machine$double.eps)
)

# S4 method for class 'markovchain'
closestReversible(
  object,
  stationaryDistribution = NULL,
  tolerance = sqrt(.Machine$double.eps)
)
```

## Arguments

- object:

  A `markovchain` object representing a finite, irreducible
  discrete-time Markov chain.

- stationaryDistribution:

  Optional numeric vector giving the stationary distribution \\\pi\\ to
  make the result reversible with respect to. It must have one strictly
  positive entry per state, sum to one (it is normalized if it does
  not), and be stationary for `object`, i.e. satisfy \\\pi P = \pi\\;
  otherwise an error is raised. The default, `NULL`, uses the chain's
  own unique stationary distribution from
  [`steadyStates`](steadyStates.md).

- tolerance:

  A single finite non-negative number, used when checking that a
  supplied `stationaryDistribution` really is stationary. It is a
  numerical-noise tolerance, not a modelling one.

## Value

A named list with four elements:

- `chain`:

  The approximating `markovchain` object \\R\\, with the same states and
  the same row/column-stochastic storage convention as `object`.

- `stationaryDistribution`:

  The \\\pi\\ used, as a named numeric vector.

- `distance`:

  The distance \\\\P-R\\\_\pi\\ actually minimized (see Details).

- `frobeniusDistance`:

  The plain Frobenius distance \\\\P-R\\\_F\\, reported for convenience.
  It is *not* the quantity being minimized.

## Details

For a row-stochastic transition matrix \\P\\ with stationary
distribution \\\pi\\, define the *time reversal* \\P^\*\\ by
\$\$P^\*\_{ij} = \frac{\pi_j P\_{ji}}{\pi_i}.\$\$ \\P^\*\\ is the
transition matrix of the same chain run backwards in time, and \\P\\ is
reversible exactly when \\P = P^\*\\. The approximation returned is the
*additive reversibilization* \$\$R = \tfrac{1}{2}\left(P +
P^\*\right),\$\$ which is stochastic, non-negative, and reversible with
respect to the same \\\pi\\ (see Details).

**In what sense is this the closest chain?** Work in the space
\\\ell^2(\pi)\\ of functions on the states with inner product \\\langle
f,g\rangle\_\pi=\sum_i \pi_i f_i g_i\\. Reversible chains are exactly
the self-adjoint operators on that space, and they form a linear
subspace. The matching Hilbert-Schmidt inner product on operators is
\$\$\langle A,B\rangle\_\pi = \sum\_{i,j} \frac{\pi_i}{\pi_j}
A\_{ij}B\_{ij}, \qquad \\A\\\_\pi^2 = \sum\_{i,j} \frac{\pi_i}{\pi_j}
A\_{ij}^2,\$\$ and \\A \mapsto A^\*\\ is an isometric involution for it.
The orthogonal projection onto the fixed points of such an involution is
the average of a point and its image, so \\R=(P+P^\*)/2\\ is the closest
\\\pi\\-reversible matrix to \\P\\ in \\\\\cdot\\\_\pi\\. The stochastic
and non-negativity constraints come for free: \\P^\*\\ has row sums
\\\sum_j \pi_j P\_{ji}/\pi_i = (\pi P)\_i/\pi_i = 1\\ because \\\pi\\ is
stationary, and both \\P\\ and \\P^\*\\ are non-negative, so the
minimizer over the subspace already lies in the set of transition
matrices and no constrained optimization is needed.

**What this function does not do.** The minimization is over reversible
chains *with \\\pi\\ held fixed*, in the \\\pi\\-weighted norm above.
Two related problems are different and are not solved here:

- Minimizing the *plain* Frobenius distance \\\\P-R\\\_F\\ with \\\pi\\
  fixed. The involution \\A\mapsto A^\*\\ is not an isometry for that
  norm, so \\R\\ is generally not its minimizer; `frobeniusDistance` is
  reported only as a descriptive figure.

- Letting the stationary distribution *vary*, i.e. finding the
  reversible chain nearest to \\P\\ over all choices of \\\pi\\. That is
  a genuinely harder constrained optimization problem, studied by
  Nielsen and Weber (2015), and it needs a numerical optimizer rather
  than a closed form. If you need it, supply candidate distributions
  through `stationaryDistribution` and compare `distance` values, or use
  a dedicated implementation.

**Properties worth knowing.** \\R\\ has the same stationary distribution
\\\pi\\ as \\P\\, and it preserves the support pattern in the
symmetrized sense: \\R\_{ij}\>0\\ whenever \\P\_{ij}\>0\\ or
\\P\_{ji}\>0\\. It may therefore allow transitions the original chain
forbids, which is inherent to making a chain reversible rather than a
defect of this construction. If `object` is already reversible, \\R=P\\
and `distance` is zero (up to rounding). Fill (1991) introduces this
construction, alongside the *multiplicative* reversibilization
\\PP^\*\\, which is a different object and is not computed here.

Only irreducibility is required, not aperiodicity. Irreducibility
guarantees both a unique \\\pi\\ and \\\pi_i\>0\\ for every state, so
the division defining \\P^\*\\ is always safe.

The implementation calls [`steadyStates`](steadyStates.md) at most once
and is then \\O(n^2)\\ in time and memory for a dense \\n\\-state
transition matrix; no eigendecomposition or optimization is involved.

## References

Fill, J. A. (1991). Eigenvalue bounds on convergence to stationarity for
nonreversible Markov chains, with an application to the exclusion
process. *The Annals of Applied Probability*, 1(1).
[doi:10.1214/aoap/1177005981](https://doi.org/10.1214/aoap/1177005981)

Nielsen, A. and Weber, M. (2015). Computing the nearest reversible
Markov chain. *Numerical Linear Algebra with Applications*, 22.
[doi:10.1002/nla.1967](https://doi.org/10.1002/nla.1967)

## See also

[`is.reversible`](is.reversible.md), [`steadyStates`](steadyStates.md),
[`is.irreducible`](is.irreducible.md)

## Examples

``` r
# A directed 3-cycle is as far from reversible as a chain gets: it only
# ever moves one way round. Its closest reversible approximation is the
# undirected random walk on the same triangle.
statesNames <- c("a", "b", "c")
cycle3 <- new("markovchain", states = statesNames,
  transitionMatrix = matrix(c(0, 1, 0,
                              0, 0, 1,
                              1, 0, 0), byrow = TRUE, nrow = 3,
                            dimnames = list(statesNames, statesNames)))
is.reversible(cycle3)
#> [1] FALSE

approximation <- closestReversible(cycle3)
approximation$chain
#> Unnamed Markov chain (closest reversible) 
#>  A  3 - dimensional discrete Markov Chain defined by the following states: 
#>  a, b, c 
#>  The transition matrix  (by rows)  is defined as follows: 
#>     a   b   c
#> a 0.0 0.5 0.5
#> b 0.5 0.0 0.5
#> c 0.5 0.5 0.0
#> 
is.reversible(approximation$chain)
#> [1] TRUE
approximation$distance
#> [1] 1.224745
```
