# Closest reversible approximation of a Markov chain

Finds the Markov chain closest to a given one among the reversible
chains, either with respect to a fixed stationary distribution or, with
`stationaryDistribution = "free"`, over all stationary distributions.

## Usage

``` r
closestReversible(
  object,
  stationaryDistribution = NULL,
  tolerance = sqrt(.Machine$double.eps),
  norm = c("pi", "frobenius"),
  nstart = 5L
)

# S4 method for class 'markovchain'
closestReversible(
  object,
  stationaryDistribution = NULL,
  tolerance = sqrt(.Machine$double.eps),
  norm = c("pi", "frobenius"),
  nstart = 5L
)
```

## Arguments

- object:

  A `markovchain` object representing a finite, irreducible
  discrete-time Markov chain.

- stationaryDistribution:

  Optional numeric vector giving the stationary distribution \\\pi\\ to
  make the result reversible with respect to, with one strictly positive
  entry per state (it is normalized if it does not sum to one). With
  `norm = "pi"` it must be stationary for `object`, i.e. satisfy \\\pi P
  = \pi\\, otherwise an error is raised; with `norm = "frobenius"` any
  positive distribution is allowed, and the result is the chain nearest
  to \\P\\ among those reversible with respect to it. The default,
  `NULL`, uses the chain's own unique stationary distribution from
  [`steadyStates`](steadyStates.md). With `norm = "frobenius"` the
  string `"free"` lets the stationary distribution vary as well (see
  Details).

- tolerance:

  A single finite non-negative number, used when checking that a
  supplied `stationaryDistribution` really is stationary (`norm = "pi"`
  only). It is a numerical-noise tolerance, not a modelling one.

- norm:

  The distance to minimize: `"pi"` (the default), the \\\pi\\-weighted
  Hilbert-Schmidt norm, whose minimizer is the additive
  reversibilization above; or `"frobenius"`, the plain Frobenius norm
  \\\\P - R\\\_F\\, which needs the quadprog package.

- nstart:

  With `stationaryDistribution = "free"`, the number of starting points
  of the optimization (at least one); ignored otherwise. The first
  starts from the stationary distribution of `object` (or the uniform
  distribution if it is not unique), the second from the uniform
  distribution, the others from random distributions drawn with the
  session's random number generator.

## Value

A named list with five elements:

- `chain`:

  The approximating `markovchain` object \\R\\, with the same states and
  the same row/column-stochastic storage convention as `object`.

- `stationaryDistribution`:

  The \\\pi\\ with respect to which \\R\\ is reversible, as a named
  numeric vector: the one used or, with `"free"`, the one found.

- `distance`:

  The distance actually minimized: \\\\P-R\\\_\pi\\ for `norm = "pi"`,
  \\\\P-R\\\_F\\ for `norm = "frobenius"` (see Details).

- `frobeniusDistance`:

  The plain Frobenius distance \\\\P-R\\\_F\\, reported for convenience;
  with `norm = "pi"` it is *not* the quantity being minimized.

- `norm`:

  The norm used.

With `stationaryDistribution = "free"` a sixth element, `starts`, gives
for each starting point the Frobenius distance reached.

## Details

For a row-stochastic transition matrix \\P\\ with stationary
distribution \\\pi\\, define the *time reversal* \\P^\*\\ by
\$\$P^\*\_{ij} = \frac{\pi_j P\_{ji}}{\pi_i}.\$\$ \\P^\*\\ is the
transition matrix of the same chain run backwards in time, and \\P\\ is
reversible exactly when \\P = P^\*\\. With the default `norm = "pi"` the
approximation returned is the *additive reversibilization* \$\$R =
\tfrac{1}{2}\left(P + P^\*\right),\$\$ which is stochastic,
non-negative, and reversible with respect to the same \\\pi\\ (see
Details). With `norm = "frobenius"` it is the reversible chain nearest
to \\P\\ in the plain Frobenius norm, the problem solved by Nielsen and
Weber (2015).

**In what sense is this the closest chain?** With the default
`norm = "pi"`, work in the space \\\ell^2(\pi)\\ of functions on the
states with inner product \\\langle f,g\rangle\_\pi=\sum_i \pi_i f_i
g_i\\. Reversible chains are exactly the self-adjoint operators on that
space, and they form a linear subspace. The matching Hilbert-Schmidt
inner product on operators is \$\$\langle A,B\rangle\_\pi = \sum\_{i,j}
\frac{\pi_i}{\pi_j} A\_{ij}B\_{ij}, \qquad \\A\\\_\pi^2 = \sum\_{i,j}
\frac{\pi_i}{\pi_j} A\_{ij}^2,\$\$ and \\A \mapsto A^\*\\ is an
isometric involution for it. The orthogonal projection onto the fixed
points of such an involution is the average of a point and its image, so
\\R=(P+P^\*)/2\\ is the closest \\\pi\\-reversible matrix to \\P\\ in
\\\\\cdot\\\_\pi\\. The stochastic and non-negativity constraints come
for free: \\P^\*\\ has row sums \\\sum_j \pi_j P\_{ji}/\pi_i = (\pi
P)\_i/\pi_i = 1\\ because \\\pi\\ is stationary, and both \\P\\ and
\\P^\*\\ are non-negative, so the minimizer over the subspace already
lies in the set of transition matrices and no constrained optimization
is needed.

**The Frobenius norm with a fixed \\\pi\\.** With `norm = "frobenius"`
the result minimizes \\\\P - R\\\_F^2 = \sum\_{i,j}(P\_{ij} -
R\_{ij})^2\\ over the transition matrices that satisfy detailed balance
\\\pi_i R\_{ij} = \pi_j R\_{ji}\\. This is the problem of Nielsen and
Weber (2015). Since \\A\mapsto A^\*\\ is not an isometry of this norm,
the additive reversibilization is generally not its solution, and the
non-negativity constraints can be active. The problem is a convex
quadratic program: with the \\n(n-1)/2\\ entries \\R\_{ij}\\, \\i\<j\\,
as unknowns, \\R\_{ji} = \pi_i R\_{ij}/\pi_j\\ and \\R\_{ii} = 1 -
\sum\_{j\ne i}R\_{ij}\\, the objective is strictly convex and the
constraints are \\R\_{ij} \ge 0\\ and \\R\_{ii} \ge 0\\. It is solved
exactly by the dual active-set method of
[`quadprog::solve.QP`](https://rdrr.io/pkg/quadprog/man/solve.QP.html),
in \\O(n^6)\\ time, which is practical up to about 50 states. \\\pi\\
need not be stationary for \\P\\; the result is reversible, hence
stationary, with respect to it.

**Letting \\\pi\\ vary.** With `stationaryDistribution = "free"` the
minimization is over all reversible chains, whatever their stationary
distribution. Every such chain can be written as \\R\_{ij} =
S\_{ij}/\sum_k S\_{ik}\\ with \\S\\ symmetric and non-negative
(\\S\_{ij} = \pi_i R\_{ij}\\, the probability flows), so the problem
becomes a smooth minimization over \\S \ge 0\\, solved with
`optim(method = "L-BFGS-B")` and the analytic gradient. The set of
reversible chains is not convex, so the result is a local minimum:
`nstart` starting points are tried, and each local solution is polished
by solving the convex problem above at its own \\\pi\\, which can only
decrease the distance. Nielsen and Weber (2015) keep \\\pi\\ fixed; this
extension is specific to the package. The distance found is never larger
than the one obtained with the stationary distribution of `object`.

**What this function does not do.** With `norm = "pi"` the minimization
holds \\\pi\\ fixed: the \\\pi\\-weighted norm itself depends on
\\\pi\\, so letting it vary is only offered for the Frobenius norm.
Also, the weighted Frobenius norm of PyDTMC's
`closest_reversible(weighted = TRUE)` is not offered. PyDTMC solves the
unweighted problem with an iterative solver, whose solutions agree with
the ones here to about \\10^{-5}\\ and have a slightly larger distance.

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

# The nearest reversible chain in the Frobenius norm (needs quadprog),
# first with the stationary distribution held fixed, then free.
if (requireNamespace("quadprog", quietly = TRUE)) {
  P <- new("markovchain", states = letters[1:4], transitionMatrix = matrix(
    c(0.1, 0.6, 0.2, 0.1, 0.05, 0.15, 0.7, 0.1,
      0.3, 0.05, 0.25, 0.4, 0.5, 0.2, 0.1, 0.2), 4, byrow = TRUE,
    dimnames = list(letters[1:4], letters[1:4])))
  closestReversible(P)$frobeniusDistance
  closestReversible(P, norm = "frobenius")$distance
  set.seed(1)
  free <- closestReversible(P, "free", norm = "frobenius")
  free$distance
  free$stationaryDistribution
}
#>          a          b          c          d 
#> 0.20622910 0.29041450 0.40357449 0.09978191 
```
