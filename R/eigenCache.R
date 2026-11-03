# Session cache of eigenvalues.
#
# Several functions need the eigenvalues of the same transition matrix: slem(), spectralGap(),
# relaxationTime() and impliedTimescales() (through .nonTrivialEigenvalueModuli()),
# autoplot(type = "eigenvalues"), the automatic choice of k in aggregateStates() and the
# topological entropy. A user who calls several of them on one chain, or summary(details = TRUE),
# used to pay for the same O(n^3) decomposition every time. Here the values are computed once and
# kept for the rest of the session.
#
# How it works:
#   * the key is the matrix itself, without dimnames: a stored entry is reused only when its
#     matrix is bitwise identical to the requested one (identical(num.eq = FALSE)), so a cached
#     result is always exactly what eigen() would return now. Comparing costs O(n^2), against the
#     O(n^3) of the decomposition;
#   * the most recently used entries are kept (least recently used out first), within a memory
#     budget: getOption("markovchain.cacheSize"), in bytes, default 64 MiB. A matrix larger than
#     the budget is never stored;
#   * options(markovchain.cache = FALSE) switches the cache off: eigen() is then called every
#     time and nothing is stored.
# Nothing in the public interface changes, and results are the same with or without the cache.

# Environment holding the entries; each entry is list(matrix, values, bytes).
.eigenCache <- new.env(parent = emptyenv())
.eigenCache$entries <- list()

# TRUE unless the user set options(markovchain.cache = FALSE).
.eigenCacheEnabled <- function() {
  !isFALSE(getOption("markovchain.cache", TRUE))
}

# Memory budget of the cache in bytes (non-numeric or negative values fall back to the default).
.eigenCacheBudget <- function() {
  budget <- getOption("markovchain.cacheSize", 64 * 1024^2)
  if (length(budget) != 1L || !is.numeric(budget) || is.na(budget) || budget < 0) {
    budget <- 64 * 1024^2
  }
  budget
}

# Empties the cache (used by the tests; harmless at any time).
.clearEigenCache <- function() {
  .eigenCache$entries <- list()
  invisible(NULL)
}

# Eigenvalues of the square matrix M, as eigen(M, only.values = TRUE)$values, taken from the cache
# when the same matrix was decomposed before in this session.
.eigenvaluesCached <- function(M) {
  M <- unname(as.matrix(M))
  if (!.eigenCacheEnabled()) {
    return(eigen(M, only.values = TRUE)$values)
  }

  entries <- .eigenCache$entries
  for (i in seq_along(entries)) {
    stored <- entries[[i]]$matrix
    if (identical(dim(stored), dim(M)) && identical(stored, M, num.eq = FALSE)) {
      # hit: move the entry to the front (most recently used) and return its values
      .eigenCache$entries <- c(entries[i], entries[-i])
      return(entries[[i]]$values)
    }
  }

  values <- eigen(M, only.values = TRUE)$values

  bytes <- as.numeric(utils::object.size(M)) + as.numeric(utils::object.size(values))
  budget <- .eigenCacheBudget()
  if (bytes <= budget) {
    entries <- c(list(list(matrix = M, values = values, bytes = bytes)), entries)
    # drop the least recently used entries until the budget is respected
    used <- cumsum(vapply(entries, function(e) e$bytes, numeric(1)))
    .eigenCache$entries <- entries[used <= budget]
  }
  values
}

#' Session cache of eigenvalues
#'
#' Functions of the package that need the eigenvalues of a transition matrix
#' (\code{\link{slem}}, \code{\link{spectralGap}},
#' \code{\link{relaxationTime}}, \code{\link{impliedTimescales}},
#' \code{\link{topologicalEntropy}}, \code{\link{aggregateStates}} with
#' \code{k = NULL}, \code{autoplot(type = "eigenvalues")} and
#' \code{summary(details = TRUE)}) compute them once per matrix and keep them
#' for the rest of the R session, so that calling several of these functions
#' on the same chain does not repeat the \eqn{O(n^3)} eigen-decomposition.
#'
#' A stored result is reused only for a matrix bitwise identical to the one it
#' was computed from (state names are ignored), so results are exactly the
#' same with or without the cache. The most recently used matrices are kept
#' within a memory budget.
#'
#' @section Options:
#' \describe{
#'   \item{\code{markovchain.cache}}{Logical; \code{FALSE} switches the cache
#'     off. Default \code{TRUE}.}
#'   \item{\code{markovchain.cacheSize}}{Memory budget of the cache in bytes,
#'     default \code{64 * 1024^2} (64 MiB). A matrix larger than the budget
#'     is not stored.}
#' }
#'
#' @examples
#' mc <- new("markovchain", states = c("a", "b", "c"),
#'   transitionMatrix = matrix(c(0.5, 0.3, 0.2,
#'                               0.1, 0.8, 0.1,
#'                               0.3, 0.3, 0.4), nrow = 3, byrow = TRUE))
#' slem(mc)          # computes the eigenvalues
#' spectralGap(mc)   # reuses them
#'
#' # switch the cache off for the rest of the session
#' old <- options(markovchain.cache = FALSE)
#' spectralGap(mc)
#' options(old)
#'
#' @name markovchainCache
NULL
