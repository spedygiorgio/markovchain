#' Print a statistical test result in a consistent format
#'
#' @keywords internal
.printTestResult <- function(name, statistic, dof, p.value, method = NULL) {
  cat(sprintf("%s\n", name))
  if (!is.null(method))
    cat(sprintf("  Method             : %s\n", method))
  cat(sprintf("  Test statistic     : %10.4f\n", statistic))
  cat(sprintf("  Degrees of freedom : %10d\n", as.integer(dof)))
  cat(sprintf("  p-value            : %10.4g\n", p.value))
}
