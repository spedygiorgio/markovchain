# Implied timescales of a Markov chain

Computes the implied relaxation timescale associated with each
non-trivial eigenvalue of a finite, irreducible discrete-time Markov
chain.

## Usage

``` r
impliedTimescales(object)

# S4 method for class 'markovchain'
impliedTimescales(object)
```

## Arguments

- object:

  A `markovchain` object representing a finite, irreducible
  discrete-time Markov chain.

## Value

A named numeric vector of length \\n-1\\ (one entry per non-trivial
eigenvalue), sorted by decreasing timescale, i.e. by decreasing
eigenvalue modulus. Names are `"tau2"`, `"tau3"`, ..., matching the
usual eigenvalue indexing \\\lambda_2,\lambda_3,\ldots\\ in decreasing
modulus. For the trivial one-state chain, a length-zero named numeric
vector is returned.

## Details

For a row-stochastic transition matrix \\P\\ with eigenvalues
\\1=\lambda_1,\lambda_2,\ldots,\lambda_n\\ (\\\|\lambda_1\|\\ the unique
unit eigenvalue of an irreducible chain), the implied timescale of
\\\lambda_k\\, \\k\>1\\, is \$\$\tau_k =
-\frac{1}{\log\|\lambda_k\|}\$\$ for \\0\<\|\lambda_k\|\<1\\. Each
\\\tau_k\\ measures how many steps the mode associated with
\\\lambda_k\\ takes to decay by a factor of \\1/e\\; larger timescales
correspond to slower-decaying, more persistent modes.

Only irreducibility is required, not aperiodicity: this is the same
convention used by [`slem`](slem.md) and
[`spectralGap`](spectralGap.md), and it lets `impliedTimescales()`
document periodic and boundary cases explicitly rather than rejecting
them:

- If \\\|\lambda_k\|\\ is (numerically) exactly \\1\\ – which happens
  for non-trivial eigenvalues of periodic chains, e.g. \\\lambda=-1\\
  for a 2-cycle – the corresponding mode never decays and `tau_k = Inf`
  is returned. This is a boundary case of the formula above (as
  \\\|\lambda\|\to 1^-\\, \\\tau\to\infty\\) that is handled explicitly
  rather than by evaluating \\-1/\log(1)\\, which is numerically `-Inf`
  rather than the mathematically correct `+Inf`.

- If \\\|\lambda_k\|\\ is (numerically) exactly \\0\\, the mode decays
  immediately and `tau_k = 0` is returned. `log(0)` evaluates to `-Inf`
  in R, so this case is already handled correctly by the formula itself
  and needs no special-casing.

The term "implied timescale" follows the Markov state model literature
in molecular kinetics, where it is additionally used, across chains
estimated at increasing lag times, as a self-consistency check on the
Markov (memoryless) approximation: implied timescales that are
approximately constant across lag times support the model, while ones
that drift indicate it should be revisited (see Prinz et al. (2011)).
Building such a lag-time comparison is left to the user, since it
requires re-estimating the chain at each lag: `impliedTimescales()`
itself only evaluates a single, already-fitted `markovchain` object.

The implementation calls [`eigen()`](https://rdrr.io/r/base/eigen.html)
with `only.values = TRUE`, so it never computes eigenvectors. Its time
complexity is \\O(n^3)\\ and its memory use is \\O(n^2)\\ for a dense
\\n\\-state transition matrix. It supports both row- and
column-stochastic storage.

## References

Swope, W. C., Pitera, J. W. and Suits, F. (2004). Describing protein
folding kinetics by molecular dynamics simulations, 1: Theory. *J. Phys.
Chem. B*, 108(21), 6571-6581.

Prinz, J.-H., Wu, H., Sarich, M., Keller, B., Senne, M., Held, M.,
Chodera, J. D., Schutte, C. and Noe, F. (2011). Markov models of
molecular kinetics: Generation and validation. *Journal of Chemical
Physics*, 134(17), 174105.

## See also

[`slem`](slem.md), [`spectralGap`](spectralGap.md),
[`is.irreducible`](is.irreducible.md)

## Examples

``` r
statesNames <- c("a", "b")
mc <- new("markovchain",
  states = statesNames,
  transitionMatrix = matrix(c(0.7, 0.3, 0.1, 0.9),
    byrow = TRUE, nrow = 2,
    dimnames = list(statesNames, statesNames)))
impliedTimescales(mc)
#>     tau2 
#> 1.957615 
```
