# Methods to generate random markov chains

# Internal helper: evaluate `expr` after set.seed(seed), restoring the
# caller's random number stream afterwards, so that a `seed` argument gives
# reproducible output without side effects on the global RNG state. With
# seed = NULL, `expr` uses (and advances) the current stream as usual.
.withLocalSeed <- function(seed, expr) {
  if (is.null(seed)) {
    return(expr)
  }
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
      seed != round(seed)) {
    stop("seed must be NULL or a single whole number.")
  }
  hadSeed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (hadSeed) {
    oldSeed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  on.exit({
    if (hadSeed) {
      assign(".Random.seed", oldSeed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  set.seed(seed)
  expr
}

# Internal helper: validate the number of states and the state names of a
# generated chain; returns list(n, states).
.generatorStates <- function(n, states, minStates = 1L) {
  if (is.null(n)) {
    if (is.null(states)) {
      stop("Provide n or states.")
    }
    n <- length(states)
  }
  if (length(n) != 1L || !is.numeric(n) || !is.finite(n) || n != round(n) ||
      n < minStates) {
    stop(paste0("n must be a single whole number of at least ", minStates, "."))
  }
  n <- as.integer(n)
  if (is.null(states)) {
    states <- as.character(seq_len(n))
  } else {
    states <- as.character(states)
    if (length(states) != n || anyNA(states) || anyDuplicated(states) ||
        any(!nzchar(states))) {
      stop(paste0("states must contain ", n,
                  " distinct, non-missing, non-empty names."))
    }
  }
  list(n = n, states = states)
}

#' Random Markov chain
#'
#' Generates a Markov chain whose transition probabilities are drawn at
#' random, optionally with a given number of zero probabilities and with
#' some probabilities fixed in advance.
#'
#' @param n The number of states. It can be omitted when \code{states} is
#'   given.
#' @param states An optional character vector of \code{n} state names.
#'   Defaults to \code{as.character(1:n)}.
#' @param zeros The number of transition probabilities, among those not
#'   fixed by \code{mask}, that are set to zero. Every row whose
#'   probabilities are not all fixed keeps at least one positive free
#'   entry, which bounds \code{zeros} from above; a larger value is an error.
#' @param mask An optional \code{n x n} matrix of fixed transition
#'   probabilities: \code{NA} marks the entries to draw at random, any other
#'   value (in \eqn{[0, 1]}) is kept as it is. In each row the fixed values
#'   must not sum to more than one. A row whose fixed values sum to one gets
#'   zero in its \code{NA} entries; in any other row the free entries share
#'   the remaining probability. With \code{byrow = FALSE} the mask is read
#'   by columns, like the transition matrix.
#' @param byrow Whether the transition matrix of the result is stored by
#'   rows (the default) or by columns.
#' @param seed An optional whole number. When given, the chain is generated
#'   after \code{set.seed(seed)} and the caller's random number stream is
#'   restored afterwards, so the result is reproducible without affecting
#'   later random draws; when \code{NULL} (the default), the current stream
#'   is used, so \code{set.seed()} beforehand works as usual.
#' @param name The \code{name} slot of the result.
#'
#' @details
#' The algorithm follows \code{MarkovChain.random()} of PyDTMC. In each row
#' not completely fixed by \code{mask}, one free entry, chosen at random, is
#' reserved to be positive. The \code{zeros} zero entries are then chosen at
#' random among the remaining free entries of the whole matrix, the other
#' free entries are drawn from a uniform distribution on \eqn{(0, 1)}, and
#' the free entries of each row are rescaled so that, together with the
#' fixed ones, they sum to one. The rows are normalised uniforms, which is
#' not the uniform distribution on the simplex; use
#' \code{\link{dirichletChain}} or draw the rows yourself if that matters.
#'
#' @return A \code{markovchain} object with \code{n} states.
#'
#' @seealso \code{\link{dirichletChain}}, \code{\link{identityChain}}
#'
#' @examples
#' randomMarkovChain(4, seed = 1)
#' # 6 of the 16 transition probabilities are zero
#' sum(randomMarkovChain(4, zeros = 6, seed = 1)@transitionMatrix == 0)
#' # state "b" moves to "a" with probability 0.5; the rest is random
#' m <- matrix(NA, 3, 3)
#' m[2, 1] <- 0.5
#' randomMarkovChain(states = c("a", "b", "c"), mask = m, seed = 1)
#'
#' @export
randomMarkovChain <- function(n, states = NULL, zeros = 0L, mask = NULL,
                              byrow = TRUE, seed = NULL,
                              name = "Random Markov chain") {
  st <- .generatorStates(if (missing(n)) NULL else n, states)
  n <- st$n
  if (length(zeros) != 1L || !is.numeric(zeros) || !is.finite(zeros) ||
      zeros < 0 || zeros != round(zeros)) {
    stop("zeros must be a single non-negative whole number.")
  }
  if (length(byrow) != 1L || !is.logical(byrow) || is.na(byrow)) {
    stop("byrow must be TRUE or FALSE.")
  }
  tol <- sqrt(.Machine$double.eps)

  if (is.null(mask)) {
    M <- matrix(NA_real_, n, n)
  } else {
    if (!is.matrix(mask) || nrow(mask) != n || ncol(mask) != n ||
        !(is.numeric(mask) || all(is.na(mask)))) {
      stop(paste0("mask must be a ", n, " x ", n, " numeric matrix."))
    }
    M <- matrix(as.numeric(mask), n, n)
    fixedValues <- M[!is.na(M)]
    if (any(!is.finite(fixedValues)) || any(fixedValues < 0) ||
        any(fixedValues > 1)) {
      stop("The values fixed by mask must be probabilities in [0, 1].")
    }
    if (!byrow) {
      M <- t(M)
    }
  }

  fixedSum <- rowSums(M, na.rm = TRUE)
  if (any(fixedSum > 1 + tol)) {
    stop("The probabilities fixed by mask sum to more than one in some row.")
  }
  full <- abs(fixedSum - 1) <= tol
  M[is.na(M) & full] <- 0
  free <- is.na(M)
  if (any(!full & rowSums(free) == 0)) {
    stop("In some row the probabilities fixed by mask sum to less than one and no entry is left free.")
  }
  maxZeros <- sum(free) - sum(!full)
  if (zeros > maxZeros) {
    stop(paste0("zeros cannot exceed ", maxZeros,
                ", the number of free entries that can be zero."))
  }

  P <- .withLocalSeed(seed, {
    reserved <- integer(0)
    for (i in which(!full)) {
      cols <- which(free[i, ])
      j <- cols[sample.int(length(cols), 1L)]
      reserved <- c(reserved, (j - 1L) * n + i)
    }
    candidates <- setdiff(which(free), reserved)
    if (zeros > 0) {
      M[candidates[sample.int(length(candidates), zeros)]] <- 0
    }
    draw <- is.na(M)
    M[draw] <- stats::runif(sum(draw))
    for (i in which(!full)) {
      s <- sum(M[i, draw[i, ]])
      M[i, draw[i, ]] <- M[i, draw[i, ]] * (1 - fixedSum[i]) / s
    }
    M
  })

  if (!byrow) {
    P <- t(P)
  }
  dimnames(P) <- list(st$states, st$states)
  new("markovchain", states = st$states, transitionMatrix = P,
      byrow = byrow, name = name)
}

#' Markov chain from a Dirichlet process
#'
#' Generates a Markov chain whose rows are drawn from a truncated Dirichlet
#' process with the stick-breaking (GEM) construction, as
#' \code{MarkovChain.dirichlet_process()} of PyDTMC does.
#'
#' @param n The number of states, at least 2. It can be omitted when
#'   \code{states} is given.
#' @param diffusion The concentration parameter \eqn{\alpha > 0} of the
#'   Dirichlet process. Small values concentrate the probability of each row
#'   on its first states; large values spread it more evenly.
#' @param states An optional character vector of \code{n} state names.
#'   Defaults to \code{as.character(1:n)}.
#' @param diagonalBias An optional positive number \eqn{\beta}. When given,
#'   a draw from \eqn{\mathrm{Beta}(\beta, 1)} is added to each diagonal
#'   entry before the row is renormalised, which makes the chain more likely
#'   to stay where it is; larger values give a stronger bias.
#' @param shiftConcentration If \code{TRUE}, the columns are reversed, so
#'   that the probability concentrates on the last states instead of the
#'   first ones.
#' @param byrow Whether the transition matrix of the result is stored by
#'   rows (the default) or by columns.
#' @param seed An optional whole number, as in
#'   \code{\link{randomMarkovChain}}.
#' @param name The \code{name} slot of the result.
#'
#' @details
#' For each row, \eqn{b_1, \ldots, b_n} are independent
#' \eqn{\mathrm{Beta}(1, \alpha)} draws and the weights are
#' \deqn{w_j = b_j \prod_{k < j} (1 - b_k),}
#' normalised to sum to one (the truncation at \eqn{n} states leaves out the
#' mass \eqn{\prod_k (1 - b_k)}). PyDTMC only accepts whole values of
#' \eqn{\alpha} between 1 and \eqn{n}; any positive value is accepted here,
#' since the construction is defined for every \eqn{\alpha > 0}.
#'
#' @return A \code{markovchain} object with \code{n} states.
#'
#' @references
#' Sethuraman, J. (1994). A constructive definition of Dirichlet priors.
#' \emph{Statistica Sinica}, 4(2), 639-650.
#'
#' @seealso \code{\link{randomMarkovChain}}
#'
#' @examples
#' dirichletChain(5, diffusion = 2, seed = 1)
#' # a chain that tends to stay in its current state
#' dirichletChain(5, diffusion = 2, diagonalBias = 5, seed = 1)
#'
#' @export
dirichletChain <- function(n, diffusion, states = NULL, diagonalBias = NULL,
                           shiftConcentration = FALSE, byrow = TRUE,
                           seed = NULL, name = "Dirichlet process chain") {
  st <- .generatorStates(if (missing(n)) NULL else n, states, minStates = 2L)
  n <- st$n
  if (missing(diffusion) || length(diffusion) != 1L ||
      !is.numeric(diffusion) || !is.finite(diffusion) || diffusion <= 0) {
    stop("diffusion must be a single positive number.")
  }
  if (!is.null(diagonalBias) &&
      (length(diagonalBias) != 1L || !is.numeric(diagonalBias) ||
       !is.finite(diagonalBias) || diagonalBias <= 0)) {
    stop("diagonalBias must be NULL or a single positive number.")
  }
  if (length(shiftConcentration) != 1L || !is.logical(shiftConcentration) ||
      is.na(shiftConcentration)) {
    stop("shiftConcentration must be TRUE or FALSE.")
  }
  if (length(byrow) != 1L || !is.logical(byrow) || is.na(byrow)) {
    stop("byrow must be TRUE or FALSE.")
  }

  P <- .withLocalSeed(seed, {
    draws <- matrix(stats::rbeta(n * n, 1, diffusion), n, n, byrow = TRUE)
    W <- t(apply(draws, 1L, function(b) b * cumprod(c(1, 1 - b[-n]))))
    W <- W / rowSums(W)
    if (shiftConcentration) {
      W <- W[, n:1, drop = FALSE]
    }
    if (!is.null(diagonalBias)) {
      diag(W) <- diag(W) + stats::rbeta(n, diagonalBias, 1)
      W <- W / rowSums(W)
    }
    W
  })

  if (!byrow) {
    P <- t(P)
  }
  dimnames(P) <- list(st$states, st$states)
  new("markovchain", states = st$states, transitionMatrix = P,
      byrow = byrow, name = name)
}
