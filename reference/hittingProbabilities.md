# Hitting probabilities for markovchain

Given a markovchain object, this function calculates the probability of
ever arriving from state i to j

## Usage

``` r
hittingProbabilities(object, targets = NULL,
  solver = c("direct", "bicgstab", "doubling"), tol = 1e-13, maxIter = 200)
```

## Arguments

- object:

  the markovchain-class object

- targets:

  optional character vector of state names: only the hitting
  probabilities *towards* these states are computed. The default,
  `NULL`, means all the states, which gives the full matrix as before.
  Every target is handled independently of the others, so the work is
  proportional to the number of targets: for a large chain, asking only
  for the states of interest is much faster than computing the whole
  matrix and subsetting it. Duplicated or unknown names are an error.

- solver:

  the method used to solve the linear system \\(I - Q) h = R\\ on the
  states whose probability is neither structurally zero nor one (see
  Details). `"direct"`, the default, is an LU factorisation: it costs
  \\O(m^3)\\ once per target and is the most accurate. `"bicgstab"` is
  an unpreconditioned BiCGSTAB iteration on the sparse system: every
  iteration costs two sparse matrix-vector products instead of a dense
  \\O(m^3)\\ step, so it is the fastest choice on large sparse chains,
  at the price of a looser residual. On a breakdown of the iteration it
  restarts from the current residual, and if the breakdown persists it
  switches to `"direct"` with a warning. `"doubling"` is the doubled
  Neumann series used by versions up to 1.2, kept for reproducibility of
  earlier results; it is also the automatic fallback of `"direct"` on a
  numerically singular system.

- tol:

  relative residual at which the iterative solvers (`"bicgstab"`,
  `"doubling"`) stop. Ignored by `"direct"`, except when it falls back
  to `"doubling"`.

- maxIter:

  maximum number of iterations of the iterative solvers: the number of
  BiCGSTAB steps, or the number of squarings of the doubled Neumann
  series. A warning is raised, and the current values returned, if the
  requested `tol` is not reached within this many iterations.

## Value

a matrix of hitting probabilities. Entry `[i, j]` is the probability of
ever arriving from state `i` to state `j` (the probability of returning,
after at least one transition, on the diagonal); for a chain with
`byrow = FALSE` the matrix is transposed, as the transition matrix is.
With `targets`, only the columns (rows if `byrow = FALSE`) of the
targets are returned, in the order given, and they coincide with those
of the full matrix.

## Details

On each target the states are first split by graph reachability: a state
that cannot reach the target has probability zero, and one that can
reach the target but no closed class outside it has probability one.
Only the remaining states need the linear system that `solver` controls,
so on chains where that split already decides every state (an
irreducible chain, for instance) all three solvers do the same
negligible amount of work. The choice matters on chains with several
closed classes, i.e. on genuine absorption probabilities.

## References

R. Vélez, T. Prieto, Procesos Estocásticos, Librería UNED, 2013

H. A. van der Vorst (1992). Bi-CGSTAB: A Fast and Smoothly Converging
Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems. *SIAM
Journal on Scientific and Statistical Computing*, 13(2), 631-644.

## Author

Ignacio Cordón

## Examples

``` r
M <- markovchain:::zeros(5)
M[1,1] <- M[5,5] <- 1
M[2,1] <- M[2,3] <- 1/2
M[3,2] <- M[3,4] <- 1/2
M[4,2] <- M[4,5] <- 1/2

mc <- new("markovchain", transitionMatrix = M)
hittingProbabilities(mc)
#>     1     2     3         4   5
#> 1 1.0 0.000 0.000 0.0000000 0.0
#> 2 0.8 0.375 0.500 0.3333333 0.2
#> 3 0.6 0.750 0.375 0.6666667 0.4
#> 4 0.4 0.500 0.250 0.1666667 0.6
#> 5 0.0 0.000 0.000 0.0000000 1.0

# only the probabilities of ever reaching the first state
hittingProbabilities(mc, targets = "1")
#>     1
#> 1 1.0
#> 2 0.8
#> 3 0.6
#> 4 0.4
#> 5 0.0

# on a large sparse chain, the iterative solver avoids the dense products
hittingProbabilities(mc, targets = "1", solver = "bicgstab")
#>     1
#> 1 1.0
#> 2 0.8
#> 3 0.6
#> 4 0.4
#> 5 0.0
```
