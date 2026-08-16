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
  desc  <- packageDescription(pkgname, libname)
  packageStartupMessage('Package:  ', desc$Package, '\n',
                        'Version:  ', desc$Version, '\n', 
                        'Date:     ', desc$Date, '\n',
                        'BugReport: ', desc$BugReports, '\n')
}

# for unloading dynamic libraries

.onUnload <- function (libpath) {
  library.dynam.unload("markovchain", libpath)
}
