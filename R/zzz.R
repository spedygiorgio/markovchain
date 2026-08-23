# Author: Giorgio
###############################################################################

#' @importFrom stats rmultinom
NULL

utils::globalVariables(c(
  "group", "label_x", "label_y", "state",
  "x", "xend", "y", "yend"
))

.onLoad <- function(libname, pkgname) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    registerS3method(
      "autoplot",
      "markovchain",
      autoplot.markovchain,
      envir = asNamespace("ggplot2")
    )
  }
}

# loading the markovchain package

.onAttach <- function(libname, pkgname) {
  desc <- utils::packageDescription(pkgname, lib.loc = libname)

  # During devtools::load_all()/devtools::test(), packageDescription() can
  # return NA if the package is loaded from sources rather than from an
  # installed library. Avoid failing package attachment in that case.
  if (is.na(desc)[1]) {
    return(invisible())
  }

  packageStartupMessage('Package:  ', desc$Package, '\n',
                        'Version:  ', desc$Version, '\n',
                        'Date:     ', desc$Date, '\n',
                        'BugReport: ', desc$BugReports, '\n')
}

# for unloading dynamic libraries

.onUnload <- function (libpath) {
  library.dynam.unload("markovchain", libpath)
}
