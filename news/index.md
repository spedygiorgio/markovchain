# Changelog

## markovchain 0.9.0

CRAN release: 2022-07-01

- Fixed CI calculations

## markovchain 0.11.1

- Added explicit structural-zero support to
  [`assessStationarity()`](../reference/statisticalTests.md).
- Statistical inference functions now return R-standard `htest` objects
  and use the standard `print.htest()` output, while retaining the
  legacy `dof` field for compatibility.
- Standardized [`assessOrder()`](../reference/statisticalTests.md) on
  the same `htest` convention, including `parameter = c(df = ...)` and
  the standard test output.
- Expanded the `statistical_analysis` section with numerical test
  results and time-homogeneity examples.
- Added published literature examples for non-rejection of the Markov
  property and for empirical-versus-theoretical transition matrices.
- Corrected structural-zero regression tests and the
  [`assessStationarity()`](../reference/statisticalTests.md)
  documentation example used by package checks.
- Corrected
  [`verifyEmpiricalToTheoretical()`](../reference/statisticalTests.md)
  to use row totals when constructing expected transition counts for its
  row-wise multinomial likelihood-ratio test. This fixes a historical
  use of column totals and changes the Kullback et al. benchmark from
  6.551795 to 6.518384, without changing its inferential conclusion.
- Added explicit `absorbingStates` support to
  [`markovchainFit()`](../reference/markovchainFit.md) for MLE fits of
  sequence data with known terminal states.

## markovchain 0.11.0

- Fix for numerical issues
- Removed deprecated .Dim methods
- Add check monotonicity
- Fixed a division-by-zero vulnerability in the C++ backend of
  [`generatorToTransitionMatrix()`](../reference/generatorToTransitionMatrix.md).
  The function now safely handles absorbing states (where the diagonal
  element is exactly 0) without returning `NaN`s or crashing.
- Optimized memory allocation in the internal C++ function `clean_nas()`
  (used by [`createSequenceMatrix()`](../reference/markovchainFit.md)
  and other fitting functions). By pre-calculating the exact number of
  valid elements before array instantiation, the package now avoids
  continuous memory reallocations, resulting in faster parsing of large
  sequences containing `NA` values.

## markovchain 0.10.3

CRAN release: 2026-02-02

Handled *R_CHECK_PACKAGES_USED_IN_DEMO*

## markovchain 0.10.2

## markovchain 0.10.1
