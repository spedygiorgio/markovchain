# Inspect the 00check.log produced by R CMD check --as-cran and decide whether
# the package is publishable on CRAN.
#
#   Rscript .github/scripts/check-cran-notes.R <00check.log> <strict: true|false>
#
# * ERRORs and WARNINGs are never accepted (check-r-package already stops on
#   them; they are re-checked here so that the script can be used alone).
# * NOTEs are accepted only if they match `accepted` below. These are the
#   NOTEs that are explained in cran-comments.md and that CRAN accepts for this
#   package. Any other NOTE is reported as an annotation; when `strict` is
#   "true" (feat/fix/breaking commits, manual runs) it makes the job fail.
#
# To accept a new, justified NOTE add a pattern (a regular expression matched
# against the whole text of the check entry) and mention it in cran-comments.md.

accepted <- c(
  # compiled code (Rcpp, RcppArmadillo, RcppParallel) makes libs/ large
  "checking installed package size",
  # src/Makevars uses $(shell ...) to link RcppParallel: GNU make is declared
  # in SystemRequirements
  "GNU make is a SystemRequirements",
  # R CMD check cannot reach the time server from some runners
  "unable to verify current time",
  # shown by --as-cran on every update when the maintainer is unchanged
  "Days since last update"
)

args <- commandArgs(trailingOnly = TRUE)
log_file <- if (length(args) >= 1) args[1] else "check/markovchain.Rcheck/00check.log"
strict <- length(args) >= 2 && identical(tolower(args[2]), "true")
if (!file.exists(log_file)) stop("check log not found: ", log_file)

lines <- readLines(log_file, warn = FALSE)
starts <- grep("^\\* ", lines)
if (length(starts) == 0L) stop("no check entries found in ", log_file)
ends <- c(starts[-1L] - 1L, length(lines))
entries <- Map(function(s, e) lines[s:e], starts, ends)

status_of <- function(x) {
  m <- regmatches(x[1L], regexpr("\\.\\.\\. .*$", x[1L]))
  if (length(m) == 0L) "" else sub("^\\.\\.\\. ", "", sub("\\s*\\[[^]]*\\]\\s*", " ", m))
}
status <- trimws(vapply(entries, function(x) {
  s <- status_of(x)
  # statuses may carry a timing such as "[8s/21s] NOTE"
  sub("^.*\\b(OK|NOTE|WARNING|ERROR)\\b.*$", "\\1", s)
}, character(1)))
text <- vapply(entries, paste, character(1), collapse = "\n")
title <- vapply(entries, function(x) sub("\\s*\\.\\.\\..*$", "", sub("^\\* ", "", x[1L])), character(1))

bad <- status %in% c("WARNING", "ERROR")
notes <- status == "NOTE"
is_accepted <- vapply(text, function(t) any(vapply(accepted, grepl, logical(1), x = t)), logical(1))
unexpected <- notes & !is_accepted

emit <- function(level, idx) {
  for (i in idx) {
    # newlines and percent signs must be escaped in workflow commands
    msg <- gsub("\n", "%0A", gsub("%", "%25", text[i], fixed = TRUE), fixed = TRUE)
    cat(sprintf("::%s title=R CMD check %s: %s::%s\n", level, status[i], title[i], msg))
  }
}
emit("error", which(bad))
emit(if (strict) "error" else "warning", which(unexpected))
emit("notice", which(notes & is_accepted))

summary_file <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(summary_file)) {
  out <- c(
    "### CRAN readiness",
    sprintf("Mode: **%s**", if (strict) "strict (feat/fix/breaking or manual run): unexpected NOTEs fail the job" else "lenient: unexpected NOTEs are reported only"),
    sprintf("- errors/warnings: %d", sum(bad)),
    sprintf("- unexpected NOTEs: %d", sum(unexpected)),
    sprintf("- accepted NOTEs: %d", sum(notes & is_accepted)),
    if (any(bad | unexpected)) c("", "```", text[bad | unexpected], "```")
  )
  cat(out, file = summary_file, sep = "\n", append = TRUE)
}

cat(sprintf("errors/warnings: %d; unexpected NOTEs: %d; accepted NOTEs: %d; strict: %s\n",
            sum(bad), sum(unexpected), sum(notes & is_accepted), strict))
if (any(bad) || (strict && any(unexpected))) {
  cat("The package is NOT publishable on CRAN as it stands.\n")
  quit(status = 1L)
}
cat("No blocking problem found.\n")
