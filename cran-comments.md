## Submission

This is an update of markovchain 1.1.1 (published on 2026-09-18), version
1.2.1. Versions 1.1.2 and 1.2 were never submitted; their changes are included
here. Sorry for the short interval since the last release. This update fixes
functions that gave wrong results:

* `predictHommc()` simulated from the wrong probabilities (matrices read by row
  instead of by column).
* `assessOrder()` returned a zero statistic for factor input.
* `fitHigherOrder()` returned its starting weights.
* `markovchainFit(byrow = FALSE)` and several other functions gave wrong results
  for chains stored by column (`byrow = FALSE`).

It also adds new functions (among others `fitMTD()`, `higherOrderLogLik()`,
`higherOrderPredict()`, `higherOrderSimulate()`, `assessIndependence()`,
`selectOrder()`, `bootstrapProperties()`, `compareChains()`,
`divergenceRate()`, `fitFunction()`, `tauchenHussey()`, `plotSequence()`),
described in NEWS.md, and the suggested package 'quadprog'. No existing
function signature has changed. One default has changed, deliberately:
`fitHigherOrder()` now estimates the weights by maximum likelihood
(`method = "mle"`); the former least squares fit is still available with
`method = "lsq"`. This is stated in NEWS.md.

## Test environments

* local Ubuntu 24.04, R 4.3.3 (with Quarto 1.7.32 for the Quarto vignettes)
* GitHub Actions: ubuntu-latest (R release, R-devel and R 4.4), windows-latest
  (R release and R-devel), macos-latest (R release, arm64)
* TO DO before submitting: win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 2 notes

* checking installed package size ... NOTE
  installed size is 18.0Mb; sub-directories of 1Mb or more: R 1.0Mb, libs 16.3Mb
  (local check; the figures vary by platform, about 20Mb on others).
  The size of libs is due to the compiled C++ code (Rcpp, RcppArmadillo,
  RcppParallel), almost all of it debugging information: the shared library
  is about 1.1Mb once stripped.

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

The local check also reported issues of the check environment only: no
inconsolata font, no GhostScript or tidy, an unsettable locale, and Suggests
not installed. The PDF manual builds without errors with the times font.

## Reverse dependencies

TO DO: run revdepcheck before submitting.
