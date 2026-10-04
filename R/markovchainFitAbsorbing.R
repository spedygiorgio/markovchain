#' @rdname markovchainFit
#' @usage markovchainFit(data, method = "mle", byrow = TRUE, nboot = 10L,
#'   laplacian = 0, name = "", parallel = FALSE, confidencelevel = 0.95,
#'   confint = TRUE, hyperparam = matrix(), sanitize = FALSE,
#'   possibleStates = character(), absorbingStates = character(),
#'   progress = FALSE, num.cores = NULL)
#' @param absorbingStates Character vector of states that are known a priori to be
#'   absorbing. The corresponding rows are set to the identity row after MLE
#'   fitting when \code{byrow = TRUE}; the corresponding columns are set to the
#'   identity column when \code{byrow = FALSE}. The argument is currently
#'   supported only for \code{method = "mle"}.
#' @param num.cores Number of threads the parallel bootstrap path uses when
#'   \code{method = "bootstrap"} and \code{parallel = TRUE}. If \code{NULL}
#'   (the default) the thread count is read from
#'   \code{getOption("RcppParallel.numThreads")} /
#'   \code{getOption("Ncpus")} / \code{OMP_NUM_THREADS} /
#'   \code{RCPP_PARALLEL_NUM_THREADS}, falling back to \code{min(2, cores)}
#'   as CRAN policy requires. Ignored when \code{parallel = FALSE}.
#' @details When \code{absorbingStates} is supplied, the declared states must have
#'   no observed outgoing transitions. This allows terminal states in censored
#'   customer journeys to be represented as absorbing states without adding
#'   artificial observations.
#' @details \code{sanitize = "absorbing"} reaches the same result without
#'   naming the states: every state that has no observed outgoing transition,
#'   which is what \code{possibleStates} typically introduces, is made
#'   absorbing. Unlike \code{absorbingStates}, it works with every
#'   \code{method}, since it only replaces the uniform row that
#'   \code{sanitize = TRUE} would have produced. As with
#'   \code{absorbingStates}, only the estimate is constrained: any confidence
#'   bounds and standard errors keep the values the unconstrained fit assigned
#'   to those rows.
#' @export
markovchainFit <- function(data, method = "mle", byrow = TRUE, nboot = 10L,
                           laplacian = 0, name = "", parallel = FALSE,
                           confidencelevel = 0.95, confint = TRUE,
                           hyperparam = matrix(), sanitize = FALSE,
                           possibleStates = character(),
                           absorbingStates = character(), progress = FALSE,
                           num.cores = NULL) {
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be TRUE or FALSE")
  }
  # When the bootstrap worker will actually run in parallel, configure
  # RcppParallel on the main thread so .markovchainFitRcpp uses the agreed
  # number of threads. Keeping this out of C++ means set.seed() and the
  # thread choice are decided from R, consistently with rmarkovchain().
  if (isTRUE(parallel) && identical(method, "bootstrap")) {
    RcppParallel::setThreadOptions(.mcDesiredThreads(num.cores))
  }
  .markovchainFitWithAbsorbingStates(
    data, method, byrow, nboot, laplacian, name, parallel,
    confidencelevel, confint, hyperparam, sanitize, possibleStates,
    absorbingStates, progress
  )
}

# Wrapper adding explicit absorbing-state constraints to the MLE fit.
.markovchainFitWithAbsorbingStates <- function(data, method, byrow, nboot,
                                               laplacian, name, parallel,
                                               confidencelevel, confint,
                                               hyperparam, sanitize,
                                               possibleStates, absorbingStates,
                                               progress = FALSE) {
  if (!is.character(absorbingStates) || anyNA(absorbingStates)) {
    stop("`absorbingStates` must be a character vector without NA values")
  }
  declaredAbsorbing <- unique(absorbingStates)
  sanitizeMode <- .sanitizeMode(sanitize)

  # Nothing to constrain: the plain C++ fit handles sanitize = FALSE/TRUE.
  if (length(declaredAbsorbing) == 0L && sanitizeMode != "absorbing") {
    return(.Call(`_markovchain_markovchainFit`, data, method, byrow, nboot,
                 laplacian, name, parallel, confidencelevel, confint,
                 hyperparam, sanitizeMode == "uniform", possibleStates,
                 progress))
  }

  # Explicitly declared absorbing states remain an MLE-only feature, as
  # documented. sanitize = "absorbing" is not restricted that way: it only
  # replaces the uniform row that sanitize = TRUE would have produced, which
  # every method supports.
  if (length(declaredAbsorbing) > 0L && !identical(method, "mle")) {
    stop("`absorbingStates` is currently supported only with method = \"mle\"")
  }

  countData <- data
  if (is.data.frame(data) && !byrow) {
    countData <- t(as.matrix(data))
  } else if (is.matrix(data) && !byrow) {
    countData <- t(data)
  }

  # Include explicitly declared absorbing states so that their rows are retained
  # both during validation and in the final fitted transition matrix.
  fitPossibleStates <- unique(c(possibleStates, declaredAbsorbing))

  counts <- createSequenceMatrix(
    countData,
    toRowProbs = FALSE,
    sanitize = FALSE,
    possibleStates = fitPossibleStates
  )

  rowTotals <- rowSums(counts)
  hasOutgoing <- declaredAbsorbing[rowTotals[declaredAbsorbing] > 0]
  if (length(hasOutgoing) > 0L) {
    stop(sprintf(
      "Declared absorbing state(s) have observed outgoing transitions: %s",
      paste(hasOutgoing, collapse = ", ")
    ))
  }

  # sanitize = "absorbing" (#213): a state with no observed outgoing
  # transition, which is what `possibleStates` typically introduces, becomes
  # absorbing rather than uniform over every state. Such states are found
  # from the observed counts and then go through the same identity-row step
  # as the declared ones, so the two routes cannot disagree.
  absorbingStates <- declaredAbsorbing
  if (sanitizeMode == "absorbing") {
    absorbingStates <- unique(c(absorbingStates,
                                rownames(counts)[rowTotals == 0]))
  }

  # The identity rows are written below, so the fit itself must not also
  # spread those rows uniformly.
  fit <- .Call(`_markovchain_markovchainFit`, data, method, byrow, nboot,
               laplacian, name, parallel, confidencelevel, confint,
               hyperparam, FALSE, fitPossibleStates, progress)

  if (length(absorbingStates) == 0L) {
    return(fit)
  }

  transitionMatrix <- fit$estimate@transitionMatrix
  if (byrow) {
    transitionMatrix[absorbingStates, ] <- 0
    transitionMatrix[cbind(absorbingStates, absorbingStates)] <- 1
  } else {
    transitionMatrix[, absorbingStates] <- 0
    transitionMatrix[cbind(absorbingStates, absorbingStates)] <- 1
  }
  fit$estimate@transitionMatrix <- transitionMatrix

  fit
}
