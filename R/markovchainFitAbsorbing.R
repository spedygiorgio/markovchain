#' @rdname markovchainFit
#' @param absorbingStates Character vector of states that are known a priori to be
#'   absorbing. The corresponding rows are set to the identity row after MLE
#'   fitting. The argument is currently supported only for \code{method = "mle"}.
#' @details When \code{absorbingStates} is supplied, the declared states must have
#'   no observed outgoing transitions. This allows terminal states in censored
#'   customer journeys to be represented as absorbing states without adding
#'   artificial observations.
NULL

# Wrapper adding explicit absorbing-state constraints to the MLE fit.
.markovchainFitWithAbsorbingStates <- function(data, method, byrow, nboot,
                                               laplacian, name, parallel,
                                               confidencelevel, confint,
                                               hyperparam, sanitize,
                                               possibleStates, absorbingStates) {
  if (!is.character(absorbingStates) || anyNA(absorbingStates)) {
    stop("`absorbingStates` must be a character vector without NA values")
  }
  absorbingStates <- unique(absorbingStates)

  if (length(absorbingStates) == 0L) {
    return(.Call(`_markovchain_markovchainFit`, data, method, byrow, nboot,
                 laplacian, name, parallel, confidencelevel, confint,
                 hyperparam, sanitize, possibleStates))
  }

  if (!identical(method, "mle")) {
    stop("`absorbingStates` is currently supported only with method = \"mle\"")
  }

  countData <- data
  if (is.data.frame(data) && !byrow) {
    countData <- t(as.matrix(data))
  } else if (is.matrix(data) && !byrow) {
    countData <- t(data)
  }

  # Include explicitly declared absorbing states so that their rows are retained
  # even when they have no observed outgoing transitions. This is consistent with
  # the existing `possibleStates` mechanism for unobserved states.
  counts <- createSequenceMatrix(
    countData,
    toRowProbs = FALSE,
    sanitize = FALSE,
    possibleStates = unique(c(possibleStates, absorbingStates))
  )

  rowTotals <- rowSums(counts)
  hasOutgoing <- absorbingStates[rowTotals[absorbingStates] > 0]
  if (length(hasOutgoing) > 0L) {
    stop(sprintf(
      "Declared absorbing state(s) have observed outgoing transitions: %s",
      paste(hasOutgoing, collapse = ", ")
    ))
  }

  fit <- .Call(`_markovchain_markovchainFit`, data, method, byrow, nboot,
               laplacian, name, parallel, confidencelevel, confint,
               hyperparam, sanitize, possibleStates)

  transitionMatrix <- fit$estimate@transitionMatrix
  transitionMatrix[absorbingStates, ] <- 0
  transitionMatrix[cbind(absorbingStates, absorbingStates)] <- 1
  fit$estimate@transitionMatrix <- transitionMatrix

  fit
}

# RcppExports.R is generated and defines the original public function. This
# source file is loaded afterwards and extends that API with absorbingStates.
markovchainFit <- function(data, method = "mle", byrow = TRUE, nboot = 10L,
                           laplacian = 0, name = "", parallel = FALSE,
                           confidencelevel = 0.95, confint = TRUE,
                           hyperparam = matrix(), sanitize = FALSE,
                           possibleStates = character(),
                           absorbingStates = character()) {
  .markovchainFitWithAbsorbingStates(
    data, method, byrow, nboot, laplacian, name, parallel,
    confidencelevel, confint, hyperparam, sanitize, possibleStates,
    absorbingStates
  )
}
