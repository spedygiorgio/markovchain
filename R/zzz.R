# Author: Giorgio
###############################################################################

#' @importFrom stats rmultinom
NULL

# Names used inside ggplot2 aes() calls are not real variables; declaring them here
# silences the 'no visible binding for global variable' note of R CMD check.
utils::globalVariables(c(
  "group", "label_x", "label_y", "state",
  "x", "xend", "y", "yend", "re", "im", "step", "probability",
  "chain", "from", "to", "label"
))

# Hook run once when the namespace is loaded (before library() attaches it).
# If ggplot2 is installed we register our autoplot() method for 'markovchain'
# objects with ggplot2's generic. ggplot2 is only suggested, so it is optional.
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

# Hook run when the package is attached with library(markovchain): it prints the
# banner with name, version and bug-report address read from DESCRIPTION.
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

# Hook run when the package is unloaded: releases the compiled C++ library.
.onUnload <- function (libpath) {
  library.dynam.unload("markovchain", libpath)
}
