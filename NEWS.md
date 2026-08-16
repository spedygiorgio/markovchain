# markovchain 0.2.0

- Corrected tests

# markovchain 0.11.1

- Added explicit structural-zero support to `assessStationarity()`.
- Standardized statistical test console output and expanded the `statistical_analysis` vignette with numerical test results and time-homogeneity examples.

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

# markovchain 0.10.0

-  Uptick pandoc requirements and handling of sparse matrices

# News for version 0.9.6

-   Handling change of pandoc requirements

# News for version 0.9.5

-   Downtick R requirements

# News for version 0.9.4

-   Corrected strange characters

# News for version 0.9.3

-   Generalized application of requireNamespace(..., quietly = TRUE)
-   Other fixes to comply to newer CRAN requirements
-   Move to MIT license

# News for version 0.9.2

-   Add RcppParallel flags to PKG_LIBS

# News for version 0.9.1

## Current changes

-   Uptick Matrics requirements and modified Changelogs

## old changes

-   2022-09-23 0.9.1 Uptick Matrix reqs

-   2022-07-01 0.9.0 Bugfix a state classification error in Rcpp

-   2022-05-21 0.8.9 Removal of Matlab package dependency

-   2021-05-7 0.8.6 Fix a bug in markovchainListFit that made confusion between lists and data.frames

-   2020-12-04 0.8.5-2 Fixing unavailable software issues and language glitches

-   2020-09-21 0.8.5-1 Coping with etm unavailability

-   2020-05-21 0.8.5 Fixed a bug in markovchainListFit that made confusion between lists and data.frames

-   2020-05-04 0.8.4.1 2022-09-23 0.9.1

-   Uptick Matrix reqs 2022-09-23 0.9.1

-   Bugfix a state classification error in Rcpp 2022-07-01 0.9.0

-   Removal of Matlab package dependency 2022-05-21 0.8.9

-   Fixing unavailable software issues and language glitches 2020-12-04 0.8.5-2

-   Coping with etm unavailability 2020-09-21 0.8.5-1

-   Fixed a bug in markovchainListFit that made confusion between lists and data.frames 2021-05-7 0.8.6

-   Add small changes to cope with upcoming R 4.0.0 2020-03-15 0.8.3

-   Improved Rcpp performance 2018-12-08 0.6.9.12

-   Deep restructuring of statistical tests 2017-03-16 0.6.8
