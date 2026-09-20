# 1.1.3

- Added `sensitivity()`, the closed-form sensitivity of the stationary distribution of a finite irreducible DTMC to a perturbation of one state's transition row, `S[l,j] = pi[k]*(Z[l,j]-pi[j])` in terms of the fundamental matrix `Z` already used by `kemenyConstant()`. Unlike a positional, entry-by-entry sensitivity, `S` is built so that `d %*% S` gives the correct first-order change in `pi` for *any* zero-sum perturbation direction `d` of row `k` (Schweitzer, 1968; Meyer, 1980; Cho & Meyer, 2001), and the common two-state case (increase one transition, decrease another to compensate) is simply the difference of two rows of `S`.
- Added `mergeWith()`, returning the convex combination `(1-gamma)*P1 + gamma*P2` of two chains defined on the same state space. States are matched by *name*, not by matrix position: `other`'s transition matrix is realigned to `object`'s state order (and converted from column-stochastic if needed) before combining, and chains on different state-name sets are rejected with a clear error rather than silently combined positionally.
- Added `closestReversible()` (closes #253), returning the reversible chain closest to a given one for a fixed stationary distribution, together with the distance it minimizes. It is the additive reversibilization `(P + P*)/2` of Fill (1991), which is the orthogonal projection of `P` onto the reversible chains in the pi-weighted Hilbert-Schmidt norm of `l^2(pi)` and needs no numerical optimization: the stochasticity and non-negativity constraints hold automatically. The documentation states explicitly what is *not* solved -- the plain-Frobenius version of the same problem, and the harder problem of letting the stationary distribution vary (Nielsen & Weber, 2015).
- `markovchainFit()` now accepts `method = "laplace"` for a list of sequences (closes #165), pooling the transition counts over the sequences and smoothing them exactly as it does for a single sequence. `method = "bootstrap"` on a list is still refused -- resampling whole sequences and resampling transitions within them are different procedures with different standard errors, and the package will not pick one silently -- but the old `"method not available for a list"` is replaced by a message that says which methods to use instead.
- Fixed `markovchainFit(..., byrow = FALSE)` returning an **invalid** `markovchain` object from every sequence/list fitting path: `mle`, `laplace` and `map` transposed the transition matrix but left the `byrow` slot at its `TRUE` default, while `bootstrap` did the mirror image (slot set to `FALSE`, matrix left row-stochastic). Such objects failed `validObject()` and were silently misread by every method that consults the slot (`steadyStates()`, `is.irreducible()`, and so on), so fits made with `byrow = FALSE` could yield wrong results downstream. The accompanying `standardError` and confidence-interval matrices now follow the estimate's orientation as well. Matrix/data.frame input is unchanged: there `byrow` describes the layout of the observations, and the fitted chain stays row-stochastic.
- Added `is.reversible()`, checking detailed balance with respect to the stationary distribution of a finite irreducible DTMC. Only irreducibility is required, not aperiodicity, since reversibility and periodicity are independent properties.
- Added `mixingTime()`, the total-variation mixing time of a finite irreducible, *aperiodic* DTMC. Unlike `slem()`/`spectralGap()`, this genuinely requires aperiodicity: a periodic chain's distribution never converges, so its mixing time is not defined and calling it on one raises a clear error instead of looping or returning a misleading value.
- Added `lazyChain()`, building the lazy chain `alpha*I + (1-alpha)*P`; a standard device (Levin & Peres, 2017) for removing periodicity, e.g. to make `mixingTime()` applicable.
- Added `subchain()` (closes #254), restricting a chain to a subset of states either as the raw principal submatrix (`method = "submatrix"`, e.g. the `Q` block used by `fundamentalMatrix()`) or as a properly renormalized Markov chain conditioned on staying inside the subset (`method = "renormalize"`); the two are explicitly distinguished rather than conflated.


# 1.1.2

- Added `slem()` and `spectralGap()` for finite irreducible DTMCs. Only irreducibility is required (not aperiodicity): a periodic chain correctly returns `SLEM = 1` and spectral gap `0` instead of being rejected, since it genuinely does not contract towards its stationary distribution.
- Added `impliedTimescales()`, returning the relaxation timescale of every non-trivial eigenvalue of a finite irreducible DTMC, with `Inf`/`0` handled explicitly at the `|lambda| = 1`/`0` boundaries.


# 1.1.1

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