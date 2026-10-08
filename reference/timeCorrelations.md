# Time correlations and time relaxations of observed sequences

`timeCorrelations` computes the time autocorrelation of an observed
sequence of states, or the time cross-correlation of two sequences, at
stationarity. `timeRelaxations` computes how the expected value of the
observable defined by a sequence evolves from a given initial
distribution. They correspond to `time_correlations()` and
`time_relaxations()` of PyDTMC.

## Usage

``` r
timeCorrelations(object, sequence1, sequence2 = NULL, timePoints = 1)

# S4 method for class 'markovchain'
timeCorrelations(object, sequence1, sequence2 = NULL, timePoints = 1)

timeRelaxations(object, sequence, initial = NULL, timePoints = 1)

# S4 method for class 'markovchain'
timeRelaxations(object, sequence, initial = NULL, timePoints = 1)
```

## Arguments

- object:

  A `markovchain` object.

- sequence1, sequence:

  A sequence of states of the chain (character vector or factor).

- sequence2:

  An optional second sequence of states. If `NULL` (the default),
  `sequence1` is used, which gives the autocorrelation.

- timePoints:

  A vector of non-negative whole numbers, the lags at which the
  quantities are computed.

- initial:

  The initial distribution: `NULL` (uniform, the default), a single
  state, or a numeric probability vector, as in
  [`redistribute`](redistribute.md).

## Value

A numeric vector with one value per element of `timePoints`, named after
them.

## Details

A sequence defines an observable \\f\\ on the states: \\f_j\\ is the
number of times state \\j\\ occurs in it. With \\f\\ from `sequence1`,
\\g\\ from `sequence2`, transition matrix \\P\\ and stationary
distribution \\\pi\\, \$\$\mathrm{timeCorrelations}(t) = \sum_i \pi_i
f_i (P^t g)\_i = E\_\pi\[f(X_0) g(X_t)\],\$\$ and, with initial
distribution \\\mu\\, \$\$\mathrm{timeRelaxations}(t) = \mu P^t f =
E\_\mu\[f(X_t)\].\$\$ For an ergodic chain, both converge as \\t\\
grows, to \\E\_\pi\[f\] E\_\pi\[g\]\\ and to \\E\_\pi\[f\]\\
respectively, at a speed governed by the second largest eigenvalue
modulus ([`slem`](slem.md)).

The powers of \\P\\ are applied by repeated multiplication, and by
repeated squaring for long lags, never through an eigendecomposition.
PyDTMC 9.0.0 switches to an eigendecomposition as soon as a lag exceeds
the number of states; since its left and right eigenvectors are not
biorthonormal when \\P\\ has complex eigenvalues, it then returns wrong
values at every lag for such chains (the tests of this function include
one, whose correct values were checked with `numpy`).

`timeCorrelations` needs a unique stationary distribution, i.e. exactly
one recurrent class, and stops otherwise (PyDTMC returns `None`).
`timeRelaxations` is defined for every chain; unlike PyDTMC, it does not
require a unique stationary distribution.

## References

Noe, F., Doose, S., Daidone, I., Loellmann, M., Sauer, M., Chodera, J.
D. and Smith, J. C. (2011). Dynamical fingerprints for probing
individual relaxation processes in biomolecular dynamics with
simulations and kinetic experiments. *Proceedings of the National
Academy of Sciences*, 108(12), 4822-4827.

## See also

[`redistribute`](redistribute.md), [`slem`](slem.md),
[`relaxationTime`](relaxationTime.md)

## Examples

``` r
statesNames <- c("a", "b", "c")
mc <- new("markovchain", states = statesNames,
  transitionMatrix = matrix(c(0.5, 0.5, 0, 0.2, 0.3, 0.5, 0.1, 0.1, 0.8),
    byrow = TRUE, nrow = 3, dimnames = list(statesNames, statesNames)))
x <- c("a", "b", "c", "c", "c", "a")
timeCorrelations(mc, x, timePoints = 0:5)
#>        0        1        2        3        4        5 
#> 6.159091 5.715909 5.594318 5.539091 5.510818 5.496064 
timeRelaxations(mc, x, initial = "a", timePoints = c(0, 1, 10, 100))
#>        0        1       10      100 
#> 2.000000 1.500000 2.338087 2.340909 
```
