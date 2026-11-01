# Markov chain approximating a continuous transition kernel by quadrature

Builds a finite-state Markov chain from a transition density (or any
non-negative transition kernel) of a continuous-state process, by
evaluating the kernel at the nodes of a quadrature rule (the Nystrom
method). It is the counterpart of PyDTMC's `fit_function`.

## Usage

``` r
fitFunction(
  f,
  size,
  quadrature = c("gauss-legendre", "gauss-chebyshev", "fejer", "gauss-hermite",
    "trapezoid", "simpson", "midpoint"),
  interval = c(0, 1),
  mean = 0,
  sd = 1,
  states = NULL,
  name = NULL
)
```

## Arguments

- f:

  A function of two arguments, `f(x, y)`, returning the (unnormalized)
  density of moving from `x` to `y`. It should be vectorized: it is
  called once with two vectors of length `size^2`. A function that is
  not vectorized also works (it is then called once per pair), more
  slowly. Values must be finite and non-negative.

- size:

  A single integer of at least `2`: the number of nodes, i.e. of states
  of the chain.

- quadrature:

  The quadrature rule, one of `"gauss-legendre"` (the default),
  `"gauss-chebyshev"`, `"fejer"`, `"gauss-hermite"`, `"trapezoid"`,
  `"simpson"` and `"midpoint"`; see Details.

- interval:

  A numeric vector `c(a, b)` with `a < b`, the interval on which the
  process lives, for every rule except `"gauss-hermite"`.

- mean, sd:

  Mean and standard deviation of the normal distribution that defines
  the nodes of `"gauss-hermite"`, which covers the whole real line;
  ignored by the other rules.

- states:

  Optional character vector of `size` distinct state names. By default
  the nodes, formatted to 4 significant digits, are used (or `"1"`, ...,
  `"size"` if that formatting is not unique).

- name:

  Optional name of the chain.

## Value

A list with three elements:

- `chain`:

  The `markovchain` object (row-stochastic).

- `states`:

  The nodes, a numeric vector named by the states, in the order of the
  states of `chain`.

- `weights`:

  The quadrature weights for the plain integral \\\int g(y)\\dy \approx
  \sum_j w_j g(y_j)\\, named by the states.

## Details

With nodes \\y_1, \ldots, y_n\\ and weights \\w_1, \ldots, w_n\\ of a
rule for \\\int g(y)\\dy\\, the transition matrix is \$\$P\_{ij} =
\frac{w_j f(y_i, y_j)}{\sum\_{k} w_k f(y_i, y_k)}.\$\$ If \\f(x,
\cdot)\\ is the transition density of the continuous process, \\\sum_j
P\_{ij} h(y_j)\\ approximates \\E\[h(Y\_{t+1}) \mid Y_t = y_i\]\\ for
every smooth \\h\\, and the normalization corrects the mass lost outside
the nodes. The kernel does not need to integrate to one: only its values
at the nodes matter.

The rules, all with positive weights:

- `"gauss-legendre"`: Gauss-Legendre nodes, exact for polynomials of
  degree \\2n - 1\\. A good default for smooth kernels on a bounded
  interval.

- `"gauss-chebyshev"`: Gauss-Chebyshev nodes of the first kind, \\t_k =
  -\cos((2k - 1)\pi / (2n))\\ on \\\[-1, 1\]\\, with the weights
  \\(\pi/n)\sqrt{1 - t_k^2}\\ of a plain integral. It is exact for
  \\g(t) = q(t)/\sqrt{1 - t^2}\\ with \\q\\ a polynomial of degree
  \\2n - 1\\, so it suits kernels with that kind of behaviour at the end
  points.

- `"fejer"`: Fejer's first rule, the same Chebyshev nodes with the
  weights that integrate exactly every polynomial of degree \\n - 1\\.
  This is the rule that PyDTMC calls `"gauss-chebyshev"`.

- `"gauss-hermite"`: Gauss-Hermite nodes for the normal distribution
  \\N(\code{mean}, \code{sd}^2)\\, for processes on the real line. The
  weights of the plain integral are the probability weights divided by
  the normal density at the nodes. With an AR(1) density this is the
  method of Tauchen and Hussey (1991); see
  [`tauchenHussey`](tauchenHussey.md).

- `"trapezoid"` and `"simpson"`: equally spaced nodes including the end
  points; Simpson's rule needs an odd `size`.

- `"midpoint"`: the centres of `size` equal cells.

Newton-Cotes rules of high order are not offered, because their weights
become negative and would give negative probabilities.

A row whose kernel values are all zero at the nodes (the process cannot
reach any node from there) stops the function with an error: PyDTMC
replaces such a row with a uniform distribution, which invents
transitions; here the interval, the number of nodes or the kernel should
be changed instead.

For the same kernel, interval and number of nodes, the
`"gauss-legendre"`, `"fejer"`, `"trapezoid"` and `"simpson"` rules give
the same matrices as PyDTMC 10.0.0 (`"gauss-legendre"`,
`"gauss-chebyshev"`, `"trapezoid-rule"` and `"simpson-rule"`); PyDTMC's
midpoint rule (`"newton-cotes"`) ignores the lower end of the interval.

## References

Tauchen, G. and Hussey, R. (1991). Quadrature-based methods for
obtaining approximate solutions to nonlinear asset pricing models.
*Econometrica*, 59(2), 371-396.
[doi:10.2307/2938261](https://doi.org/10.2307/2938261)

Golub, G. H. and Welsch, J. H. (1969). Calculation of Gauss quadrature
rules. *Mathematics of Computation*, 23(106), 221-230.
[doi:10.2307/2004418](https://doi.org/10.2307/2004418)

Fejer, L. (1933). Mechanische Quadraturen mit positiven Cotesschen
Zahlen. *Mathematische Zeitschrift*, 37, 287-309.
[doi:10.1007/BF01474575](https://doi.org/10.1007/BF01474575)

## See also

[`tauchenHussey`](tauchenHussey.md), [`tauchen`](tauchen.md),
[`rouwenhorst`](rouwenhorst.md)

## Examples

``` r
# A process on [0, 1] that moves towards 0.6 times its current value,
# with a Gaussian kernel truncated to the interval.
kernel <- function(x, y) dnorm(y, mean = 0.6 * x + 0.2, sd = 0.15)
out <- fitFunction(kernel, size = 7, interval = c(0, 1))
out$chain
#> Quadrature approximation (gauss-legendre, 7 nodes) 
#>  A  7 - dimensional discrete Markov Chain defined by the following states: 
#>  0.02545, 0.12920, 0.29710, 0.50000, 0.70290, 0.87080, 0.97460 
#>  The transition matrix  (by rows)  is defined as follows: 
#>              0.02545      0.12920     0.29710   0.50000     0.70290
#> 0.02545 8.360400e-02 3.412097e-01 0.473186072 0.0991878 0.002783226
#> 0.12920 4.333595e-02 2.357194e-01 0.520178020 0.1912008 0.009407870
#> 0.29710 1.090360e-02 9.437633e-02 0.441444909 0.4024013 0.049102801
#> 0.50000 1.155157e-03 1.753256e-02 0.203378028 0.5558685 0.203378028
#> 0.70290 6.413373e-05 1.706878e-03 0.049102801 0.4024013 0.441444909
#> 0.87080 3.643080e-06 1.542875e-04 0.009407870 0.1912008 0.520178020
#> 0.97460 5.081911e-07 2.868414e-05 0.002783226 0.0991878 0.473186072
#>              0.87080      0.97460
#> 0.02545 2.868414e-05 5.081911e-07
#> 0.12920 1.542875e-04 3.643080e-06
#> 0.29710 1.706878e-03 6.413373e-05
#> 0.50000 1.753256e-02 1.155157e-03
#> 0.70290 9.437633e-02 1.090360e-02
#> 0.87080 2.357194e-01 4.333595e-02
#> 0.97460 3.412097e-01 8.360400e-02
#> 
# expected next value from each node, against the continuous one
cbind(chain = as.numeric(out$chain@transitionMatrix %*% out$states),
      target = 0.6 * out$states + 0.2)
#>             chain    target
#> 0.02545 0.2383721 0.2152676
#> 0.12920 0.2884502 0.2775406
#> 0.29710 0.3808824 0.3782465
#> 0.50000 0.5000000 0.5000000
#> 0.70290 0.6191176 0.6217535
#> 0.87080 0.7115498 0.7224594
#> 0.97460 0.7616279 0.7847324
```
