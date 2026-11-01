# Discretize an AR(1) process into a Markov chain (Tauchen-Hussey method)

Approximates the stationary first-order autoregressive process \$\$y_t =
(1-\rho)\alpha + \rho y\_{t-1} + \varepsilon_t, \qquad \varepsilon_t
\overset{\mathrm{iid}}{\sim} \mathcal N(0,\sigma^2)\$\$ by a
finite-state Markov chain on Gauss-Hermite nodes, following Tauchen and
Hussey (1991), with the choice of the standard deviation of the nodes
proposed by Floden (2008) as an option.

## Usage

``` r
tauchenHussey(
  alpha,
  sigma,
  rho,
  size,
  sigmaWeight = c("conditional", "floden", "unconditional")
)
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

- sigmaWeight:

  How the standard deviation \\\hat\sigma\\ of the normal distribution
  defining the nodes is chosen: `"conditional"` (the default, Tauchen
  and Hussey's original choice) uses \\\hat\sigma = \sigma\\;
  `"unconditional"` uses \\\hat\sigma = \sigma_y =
  \sigma/\sqrt{1-\rho^2}\\; `"floden"` uses Floden's compromise
  \\\hat\sigma = w\sigma + (1-w)\sigma_y\\ with \\w = 1/2 + \rho/4\\.

## Value

A named list with two elements, `chain` and `states`, in the same form
as returned by [`tauchen`](tauchen.md).

## Details

The nodes \\y_1, \ldots, y_n\\ and probability weights \\\omega_1,
\ldots, \omega_n\\ are those of the \\n\\-point Gauss-Hermite rule for
\\N(\alpha, \hat\sigma^2)\\, and \$\$P\_{ij} \propto \omega_j
\frac{\phi(y_j;\\ (1-\rho)\alpha + \rho y_i,\\ \sigma)}{\phi(y_j;\\
\alpha,\\ \hat\sigma)},\$\$ normalized by rows, with \\\phi(\cdot; m,
s)\\ the normal density. It is [`fitFunction`](fitFunction.md) with the
AR(1) transition density and `quadrature = "gauss-hermite"`.

The method is accurate for moderate persistence, but with \\\hat\sigma =
\sigma\\ the nodes are too concentrated when \\\rho\\ is close to one
and the chain understates the unconditional variance; with \\\hat\sigma
= \sigma_y\\ it overstates the persistence and understates the
conditional variance. Floden (2008) shows that his weighted choice is a
good compromise; for very persistent processes
[`rouwenhorst`](rouwenhorst.md) is usually preferable (Kopecky and Suen
2010). The vignette reproduces the conditional and unconditional
standard deviations of Floden's Table 1.

## References

Tauchen, G. and Hussey, R. (1991). Quadrature-based methods for
obtaining approximate solutions to nonlinear asset pricing models.
*Econometrica*, 59(2), 371-396.
[doi:10.2307/2938261](https://doi.org/10.2307/2938261)

Floden, M. (2008). A note on the accuracy of Markov-chain approximations
to highly persistent AR(1) processes. *Economics Letters*, 99(3),
516-520.

Kopecky, K. A. and Suen, R. M. H. (2010). Finite state Markov-chain
approximations to highly persistent processes. *Review of Economic
Dynamics*, 13(3), 701-714.

## See also

[`fitFunction`](fitFunction.md), [`tauchen`](tauchen.md),
[`rouwenhorst`](rouwenhorst.md)

## Examples

``` r
out <- tauchenHussey(alpha = 0, sigma = sqrt(0.03), rho = 0.95, size = 5,
                     sigmaWeight = "floden")
out$chain
#> Tauchen-Hussey AR(1) Approximation (rho = 0.95) 
#>  A  5 - dimensional discrete Markov Chain defined by the following states: 
#>  -0.7809, -0.3706, 0.0000, 0.3706, 0.7809 
#>  The transition matrix  (by rows)  is defined as follows: 
#>              -0.7809      -0.3706       0.0000       0.3706       0.7809
#> -0.7809 9.207167e-01 7.920494e-02 7.839379e-05 8.692134e-10 1.546272e-17
#> -0.3706 4.762712e-02 8.485696e-01 1.036614e-01 1.418607e-04 5.226711e-10
#> 0.0000  3.971920e-05 8.734415e-02 8.252323e-01 8.734415e-02 3.971920e-05
#> 0.3706  5.226711e-10 1.418607e-04 1.036614e-01 8.485696e-01 4.762712e-02
#> 0.7809  1.546272e-17 8.692134e-10 7.839379e-05 7.920494e-02 9.207167e-01
#> 
# unconditional standard deviation of the chain (Floden 2008, Table 1: 0.4792)
pi <- as.numeric(steadyStates(out$chain))
sqrt(sum(pi * (out$states - sum(pi * out$states))^2))
#> [1] 0.4792287
```
