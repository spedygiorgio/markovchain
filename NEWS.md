# 1.1.1

- Added `slem()` and `spectralGap()` for finite irreducible DTMCs. Only irreducibility is required (not aperiodicity): a periodic chain correctly returns `SLEM = 1` and spectral gap `0` instead of being rejected, since it genuinely does not contract towards its stationary distribution.
- Added `impliedTimescales()`, returning the relaxation timescale of every non-trivial eigenvalue of a finite irreducible DTMC, with `Inf`/`0` handled explicitly at the `|lambda| = 1`/`0` boundaries.
- Added `entropyRate()` for finite irreducible DTMCs, with configurable logarithm base and support for row- and column-stochastic storage.
- Added `kemenyConstant()` for finite irreducible DTMCs, with support for both row- and column-stochastic storage and tests for the zero-diagonal hitting-time convention.
- Hardened public R and native C++ entry points against invalid indices, negative sizes, integer overflow, and unbounded simulations.


# 1.1.0

- Added exact and approximate lumpability tools.
- Made `hittingProbabilities()` robust to arbitrarily small positive transitions by combining graph-based zero/one classification with a relative-residual Neumann iteration.
- Added `fundamentalMatrix()` for finite absorbing chains, including support for column-stochastic storage.
- Improved `absorptionProbabilities()` by solving `(I - Q) B = R` directly instead of explicitly forming `(I - Q)^{-1}`.

# markovchain 0.9.0

- Fixed CI calculations

# markovchain 0.11.1

- Added explicit structural-zero support to `assessStationarity()`.
- Statistical inference functions now return R-standard `htest` objects and use the standard `print.htest()` output, while retaining the legacy `dof` field for compatibility.
- Standardized `assessOrder()` on the same `htest` convention, including `parameter = c(df = ...)` and the standard test output.
- Expanded the `statistical_analysis` section with numerical test results and time-homogeneity examples.
- Added published literature examples for non-rejection of the Markov property and for empirical-versus-theoretical transition matrices.
- Corrected structural-zero regression tests and the `assessStationarity()` documentation example used by package checks.
- Corrected `verifyEmpiricalToTheoretical()` to use row totals when constructing expected transition counts for its row-wise multinomial likelihood-ratio test. This fixes a historical use of column totals and changes the Kullback et al. benchmark from 6.551795 to 6.518384, without changing its inferential conclusion.
- Added explicit `absorbingStates` support to `markovchainFit()` for MLE fits of sequence data with known terminal states.

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