# ggplot2 uses non-standard evaluation for aesthetics. Declare the data-column
# names used by autoplot.markovchain so R CMD check does not report them as
# undefined global variables.
utils::globalVariables(c(
  "group", "label_x", "label_y", "state", "x", "xend", "y", "yend"
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
