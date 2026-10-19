library(testthat)
library(markovchain)
# TEMPORARY (diagnostic branch, not to be merged): print the failures as
# GitHub annotations, because the job logs cannot be read from the session.
res <- test_check("markovchain", stop_on_failure = FALSE)
nbad <- 0L
for (t in res) for (e in t$results) {
  if (inherits(e, c("expectation_failure", "expectation_error"))) {
    nbad <- nbad + 1L
    msg <- paste0(t$file, " | ", t$test, " | ", conditionMessage(e))
    cat("\n::error::", gsub("\r?\n", "%0A", substr(msg, 1, 3000)), "\n", sep = "")
  }
}
if (nbad > 0L) stop(nbad, " test failures")
