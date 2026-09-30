# Discretize an AR(1) process into a Markov chain (Tauchen's method)

Approximates the stationary first-order autoregressive process \$\$y_t =
(1-\rho)\alpha + \rho y\_{t-1} + \varepsilon_t, \qquad \varepsilon_t
\overset{\mathrm{iid}}{\sim} \mathcal N(0,\sigma^2)\$\$ by a
finite-state Markov chain on an evenly spaced grid, following Tauchen
(1986).

## Usage

``` r
tauchen(alpha, sigma, rho, size, k = 3)
```

## Arguments

- alpha:

  A single finite number: the unconditional mean of the process.

- sigma:

  A single finite positive number: the standard deviation of the
  innovation \\\varepsilon_t\\.

- rho:

  A single number in \\(-1,1)\\: the autocorrelation (persistence) of
  the process.

- size:

  A single integer of at least `2`: the number of grid points (states)
  of the discretized chain.

- k:

  A single positive number, the half-width of the grid in units of the
  process's unconditional standard deviation
  \\\sigma_y=\sigma/\sqrt{1-\rho^2}\\. The default, `3`, follows Tauchen
  (1986)'s own recommendation and covers the great majority of the
  stationary distribution's mass for typical `rho`.

## Value

A named list with two elements:

- `chain`:

  The discretized `markovchain` object, with state names equal to the
  grid values of \\y\\ formatted to 4 significant digits.

- `states`:

  The numeric grid of \\y\\-values themselves, in the same order as
  `chain`'s states. Returning the actual levels alongside the chain,
  rather than only generic state labels `"1"`, `"2"`, ..., is
  deliberate: the whole point of discretizing an AR(1) process is
  usually to do further numeric work with the levels (e.g. plugging them
  into a pricing formula), and re-deriving the grid from `alpha`,
  `sigma`, `rho` and `k` a second time by hand is both extra work and a
  place for an off-by-one or rounding mismatch to creep in.

## Details

The grid is \\n=\code{size}\\ evenly spaced points \\y_1\<\cdots\<y_n\\
spanning \\\[\alpha-k\sigma_y,\\ \alpha+k\sigma_y\]\\, with half-spacing
\\w=(y_n-y_1)/(2(n-1))\\. Writing \\\Phi\\ for the standard normal CDF,
the transition probabilities from grid point \\y_i\\ are \$\$P\_{i1} =
\Phi\\\left(\frac{y_1-(1-\rho)\alpha-\rho y_i+w}{\sigma}\right),\$\$
\$\$P\_{in} = 1-\Phi\\\left(\frac{y_n-(1-\rho)\alpha-\rho
y_i-w}{\sigma}\right),\$\$ \$\$P\_{ij} =
\Phi\\\left(\frac{y_j-(1-\rho)\alpha-\rho y_i+w}{\sigma}\right) -
\Phi\\\left(\frac{y_j-(1-\rho)\alpha-\rho y_i-w}{\sigma}\right), \quad
1\<j\<n,\$\$ i.e. the probability that \\y_t\\ (a normal draw centred at
the AR(1) conditional mean) lands in the half-open bin around \\y_j\\,
with the two end bins extended to \\\pm\infty\\ so that rows sum to
exactly \\1\\.

Tauchen's method is simple and fast (\\O(n^2)\\ normal CDF evaluations)
but the grid width is fixed by `k` regardless of `size`: for a coarse
grid (small `size`) it under-resolves the bulk of the distribution, and
for `rho` close to \\\pm1\\ the true unconditional variance is large and
sensitive to `k`. See [`rouwenhorst`](rouwenhorst.md) for an alternative
that tends to match the persistence of near-unit-root processes more
accurately and needs no arbitrary grid-width parameter.

## References

Tauchen, G. (1986). Finite state markov-chain approximations to
univariate and vector autoregressions. *Economics Letters*, 20(2),
177-181.

## See also

[`rouwenhorst`](rouwenhorst.md)

## Examples

``` r
out <- tauchen(alpha = 0, sigma = 1, rho = 0.9, size = 5)
out$states
#>    -6.882    -3.441     0.000     3.441     6.882 
#> -6.882472 -3.441236  0.000000  3.441236  6.882472 
out$chain
#> Tauchen AR(1) Approximation (rho = 0.9) 
#>  A  5 - dimensional discrete Markov Chain defined by the following states: 
#>  -6.882, -3.441, 0.000, 3.441, 6.882 
#>  The transition matrix  (by rows)  is defined as follows: 
#>              -6.882       -3.441        0.000        3.441        6.882
#> -6.882 8.490508e-01 1.509454e-01 3.845556e-06 1.221245e-15 0.000000e+00
#> -3.441 1.947373e-02 8.961920e-01 8.433358e-02 7.260019e-07 1.110223e-16
#> 0.000  1.222580e-07 4.265996e-02 9.146798e-01 4.265996e-02 1.222580e-07
#> 3.441  7.346963e-17 7.260019e-07 8.433358e-02 8.961920e-01 1.947373e-02
#> 6.882  3.459031e-30 1.237828e-15 3.845556e-06 1.509454e-01 8.490508e-01
#> 
# The chain's own stationary variance should be close to the AR(1)'s
# theoretical unconditional variance sigma^2 / (1 - rho^2).
pi <- as.numeric(steadyStates(out$chain))
sum(pi * (out$states - sum(pi * out$states))^2)
#> [1] 8.478635
1 / (1 - 0.9^2)
#> [1] 5.263158
```
