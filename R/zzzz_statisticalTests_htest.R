# Standardize statistical inference results to R's htest convention.

.asHtest <- function(result, method, data.name) {
  statistic.name <- if (identical(method, "Pearson")) "X-squared" else "G-squared"
  result$statistic <- setNames(as.numeric(result$statistic), statistic.name)
  result$parameter <- setNames(as.numeric(result$dof), "df")
  result$dof <- NULL
  result$method <- method
  result$data.name <- data.name
  class(result) <- c("htest", setdiff(class(result), "htest"))
  result
}

# Preserve the implementations defined earlier in the package load order.
.verifyMarkovProperty <- verifyMarkovProperty
.verifyEmpiricalToTheoretical <- verifyEmpiricalToTheoretical
.verifyHomogeneity <- verifyHomogeneity
.assessStationarity <- assessStationarity

verifyMarkovProperty <- function(sequence, method = c("G", "Pearson", "simulation"),
                                 B = 9999, seed = NULL, verbose = TRUE) {
  method <- match.arg(method)
  result <- .verifyMarkovProperty(sequence, method = method, B = B,
                                  seed = seed, verbose = FALSE)
  result <- .asHtest(result, method,
                     deparse(substitute(sequence)))
  if (verbose) print(result)
  invisible(result)
}

verifyEmpiricalToTheoretical <- function(data, object,
                                         method = c("G", "Pearson", "simulation"),
                                         B = 9999, seed = NULL, verbose = TRUE) {
  method <- match.arg(method)
  result <- .verifyEmpiricalToTheoretical(data, object, method = method,
                                          B = B, seed = seed, verbose = FALSE)
  result <- .asHtest(result, method,
                     deparse(substitute(data)))
  if (verbose) print(result)
  invisible(result)
}

verifyHomogeneity <- function(inputList,
                              method = c("G", "Pearson", "simulation"),
                              B = 9999, seed = NULL, verbose = TRUE) {
  method <- match.arg(method)
  result <- .verifyHomogeneity(inputList, method = method, B = B,
                               seed = seed, verbose = FALSE)
  result <- .asHtest(result, method,
                     deparse(substitute(inputList)))
  if (verbose) print(result)
  invisible(result)
}

assessStationarity <- function(sequence, nblocks, structural.zeros = NULL,
                               verbose = TRUE) {
  result <- .assessStationarity(
    sequence, nblocks = nblocks, structural.zeros = structural.zeros,
    verbose = FALSE
  )
  result <- .asHtest(
    result, "Pearson",
    deparse(substitute(sequence))
  )
  result$method <- "Pearson's Chi-squared test for time-homogeneity"
  if (verbose) print(result)
  invisible(result)
}
