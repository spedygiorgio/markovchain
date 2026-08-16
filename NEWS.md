# markovchain 0.2.0

- Corrected tests

# markovchain 0.11.1

- Added explicit structural-zero support to `assessStationarity()`.
- Statistical inference functions now return R-standard `htest` objects and use the standard `print.htest()` output, while retaining the legacy `dof` field for compatibility.
- Standardized `assessOrder()` on the same `htest` convention, including `parameter = c(df = ...)` and the standard test output.
- Expanded the `statistical_analysis` vignette with numerical test results and time-homogeneity examples.
- Added published real-data examples from Collins (1974) and Tran and Carmichael (2012), including reconstruction of the latter's steel-fabrication transition-flow data and comparison with the package homogeneity test.
- Corrected structural-zero regression tests and the `assessStationarity()` documentation example used by package checks.

# markovchain 0.11.0

- Fix for numerical issues
- Removed deprecated .Dim methods
- Add check monotonicity
- Fixed a division-by-zero vulnerability in the C++ backend of `generatorToTransitionMatrix()`. The function now safely handles absorbing states (where the diagonal element is exactly 0) without returning `NaN`s or crashing.
- Optimized memory allocation in the internal C++ function `clean_nas()` (used by `createSequenceMatrix()` and other fitting functions). By pre-calculating the exact number of valid elements before array instantiation, the package now avoids continuous memory reallocations, resulting in faster parsing of large sequences containing `NA` values.

# markovchain 0.10.3

Handled _R_CHECK_PACKAGES_USED_IN_DEMO_

# markovchain 0.10.2

# markovchain 0.10.1