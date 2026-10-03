## Submission

This is an update of markovchain 1.1.1 (published on 2026-09-18). Version 1.1.2
was never submitted; its changes are included here. Sorry for the short interval
since the last release. This update fixes functions that gave wrong results:

* `predictHommc()` simulated from the wrong probabilities (matrices read by row
  instead of by column).
* `assessOrder()` returned a zero statistic for factor input.
* `fitHigherOrder()` returned its starting weights.
* `markovchainFit(byrow = FALSE)` returned invalid objects.

It also adds new functions (`fitMTD()`, `higherOrderLogLik()`,
`higherOrderPredict()`, `higherOrderSimulate()`, `assessIndependence()`, among
others), described in NEWS.md. No existing function signature has changed.

## Test environments

* local Ubuntu 24.04, R 4.3.3 (with Quarto 1.7.32 for the Quarto vignettes)
* GitHub Actions: ubuntu-latest (R release and R-devel), windows-latest (R release)
* TO DO before submitting: win-builder (R-devel) and R-hub / macOS builder

## R CMD check results

0 errors | 0 warnings | 2 notes

* checking installed package size ... NOTE
  installed size is 15.5Mb; sub-directories of 1Mb or more: doc 1.3Mb, libs 12.8Mb.
  The size of libs is due to the compiled C++ code (Rcpp, RcppArmadillo, RcppParallel).

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

The local check also reported issues of the check environment only: no
inconsolata font, no GhostScript or tidy, an unsettable locale, and Suggests
not installed. The PDF manual builds without errors with the times font.

## Reverse dependencies

TO DO: run revdepcheck before submitting.
