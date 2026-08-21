#' Test whether a partition is lumpable for a Markov chain
#'
#' @param object a `markovchain` object.
#' @param partition a named list defining the macro-states of the partition.
#' @return A logical value indicating whether the partition is exactly lumpable.
#' @export
setGeneric("is.lumpable", function(object, partition) standardGeneric("is.lumpable"))

#' Lump a Markov chain according to a partition
#'
#' @param object a `markovchain` object.
#' @param partition a named list defining the macro-states of the partition.
#' @param force logical; if `TRUE`, perform weighted approximate lumping when
#'   the partition is not exactly lumpable.
#' @return A lumped `markovchain` object.
#' @export
setGeneric("lump", function(object, partition, force = FALSE) standardGeneric("lump"))

#' Automatically construct a lumped Markov chain
#'
#' @param object a `markovchain` object.
#' @param k the requested number of macro-states.
#' @return A list containing the generated partition and lumped chain.
#' @export
setGeneric("autoLump", function(object, k) standardGeneric("autoLump"))

# Lumpability methods are defined here so that they are loaded after
# probabilistic.R. The helpers normalize column-stochastic chains before
# passing matrices to the row-stochastic lumpability kernels.

.get_partition_indices <- function(state_names, partition) {
  if (!is.list(partition) || length(partition) < 2L) {
    stop("partition must be a list containing at least two macro-states.")
  }
  if (is.null(names(partition)) || anyNA(names(partition)) ||
      any(!nzchar(names(partition))) || anyDuplicated(names(partition))) {
    stop("partition must be a named list with unique, non-empty macro-state names.")
  }

  part_idx <- lapply(partition, function(x) {
    if (!is.character(x) || length(x) == 0L) {
      stop("Each macro-state in partition must contain at least one state name.")
    }
    idx <- match(x, state_names)
    if (anyNA(idx)) {
      stop("Invalid partition: Some states in the partition do not exist in the Markov chain.")
    }
    idx - 1L
  })

  all_idx <- unlist(part_idx, use.names = FALSE)
  if (length(all_idx) != length(state_names) ||
      length(unique(all_idx)) != length(state_names)) {
    stop("Invalid partition: The partition must contain all states of the Markov chain exactly once (no duplicates, no omissions).")
  }

  part_idx
}

.lumpability_row_matrix <- function(object) {
  P <- object@transitionMatrix
  if (!object@byrow) {
    P <- t(P)
  }
  P
}

setMethod("is.lumpable", "markovchain", function(object, partition) {
  part_idx <- .get_partition_indices(states(object), partition)
  .is_lumpable_cpp(.lumpability_row_matrix(object), part_idx)
})

setMethod("lump", "markovchain", function(object, partition, force = FALSE) {
  part_idx <- .get_partition_indices(states(object), partition)
  P <- .lumpability_row_matrix(object)

  if (!is.logical(force) || length(force) != 1L || is.na(force)) {
    stop("force must be a single non-missing logical value.")
  }

  if (!force && !.is_lumpable_cpp(P, part_idx)) {
    stop("The Markov chain is not exactly lumpable. Use force = TRUE to perform an approximate weighted lumping.")
  }

  st <- steadyStates(object)
  w <- if (nrow(st) > 0L) {
    as.numeric(st[1, ])
  } else {
    rep(1 / ncol(P), ncol(P))
  }

  P_lumped <- .lump_cpp(P, part_idx, w)
  rownames(P_lumped) <- names(partition)
  colnames(P_lumped) <- names(partition)

  new(
    "markovchain",
    states = names(partition),
    transitionMatrix = P_lumped,
    byrow = TRUE,
    name = paste(object@name, "(Lumped)")
  )
})

setMethod("autoLump", "markovchain", function(object, k) {
  P <- .lumpability_row_matrix(object)
  state_names <- states(object)

  if (length(k) != 1L || !is.numeric(k) || is.na(k) ||
      k != as.integer(k) || k < 2L || k >= nrow(P)) {
    stop("The number of macro-states 'k' must be an integer between 2 and the number of states - 1.")
  }
  k <- as.integer(k)

  eig <- eigen(P)
  V_eig <- Re(eig$vectors[, seq_len(k), drop = FALSE])

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  set.seed(42)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  clust <- kmeans(V_eig, centers = k, nstart = 25)

  partition <- vector("list", k)
  for (i in seq_len(k)) {
    partition[[paste0("Macro_", i)]] <- state_names[clust$cluster == i]
  }

  list(
    partition = partition,
    lumped_chain = lump(object, partition, force = TRUE)
  )
})
