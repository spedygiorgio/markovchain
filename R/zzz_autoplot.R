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
