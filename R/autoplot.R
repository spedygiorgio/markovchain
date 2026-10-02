#' Plot a Markov chain with ggplot2
#'
#' Creates a ggplot2 representation of a discrete-time Markov chain. Nodes are
#' states and directed edges represent transitions with positive probability.
#' Communicating classes are shown with different node fills.
#'
#' @param object An object of class `markovchain`.
#' @param threshold Minimum transition probability to display. Defaults to 0.
#' @param show_probabilities Logical; whether to label edges with transition
#'   probabilities. Defaults to `TRUE`.
#' @param digits Number of digits used to format transition probabilities.
#' @param node_size Size of state nodes in the ggplot2 plot.
#' @param edge_width Minimum width multiplier for transition edges.
#' @param type Type of plot: \code{"graph"} (default) draws the transition
#'   graph; \code{"eigenvalues"} draws the eigenvalues of the transition matrix
#'   in the complex plane together with the unit circle; \code{"flow"} draws
#'   the evolution of the distribution over the states, as computed by
#'   \code{\link{redistribute}}; \code{"comparison"} compares \code{object}
#'   with the chains in \code{other}.
#' @param steps Number of steps of the \code{"flow"} plot. Defaults to 20.
#'   Ignored for the other types.
#' @param initial Initial distribution of the \code{"flow"} plot, see
#'   \code{\link{redistribute}}. Defaults to the uniform distribution.
#'   Ignored for the other types.
#' @param other For \code{type = "comparison"}: a \code{markovchain} object
#'   or a (possibly named) list of them, to be compared with \code{object}.
#'   All chains must be defined on the same set of states, which are matched
#'   by name, so their order may differ. Ignored for the other types.
#' @param what For \code{type = "comparison"}: \code{"transition"} (default)
#'   draws the transition matrices side by side as heatmaps on a common
#'   probability scale; \code{"stationary"} compares the stationary
#'   distributions as grouped bars (every chain must then be irreducible).
#'   Ignored for the other types.
#' @param ... Currently unused, reserved for future extensions.
#'
#' @return A ggplot object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   weather <- matrix(c(0.7, 0.2, 0.1,
#'                       0.3, 0.4, 0.3,
#'                       0.2, 0.45, 0.35),
#'                     nrow = 3, byrow = TRUE,
#'                     dimnames = list(c("sunny", "cloudy", "rain"),
#'                                     c("sunny", "cloudy", "rain")))
#'   mc <- new("markovchain", states = rownames(weather),
#'             transitionMatrix = weather, name = "Weather")
#'   ggplot2::autoplot(mc)
#'   ggplot2::autoplot(mc, type = "eigenvalues")
#'   ggplot2::autoplot(mc, type = "flow", steps = 10, initial = "rain")
#'
#'   # compare with a "stickier" version of the same chain
#'   sticky <- lazyChain(mc, alpha = 0.5)
#'   sticky@name <- "Lazy weather"
#'   ggplot2::autoplot(mc, type = "comparison", other = sticky)
#'   ggplot2::autoplot(mc, type = "comparison", other = sticky,
#'                     what = "stationary")
#' }
autoplot.markovchain <- function(object,
                                  threshold = 0,
                                  show_probabilities = TRUE,
                                  digits = 2,
                                  node_size = 6,
                                  edge_width = 1,
                                  type = c("graph", "eigenvalues", "flow",
                                           "comparison"),
                                  steps = 20,
                                  initial = NULL,
                                  other = NULL,
                                  what = c("transition", "stationary"),
                                  ...) {
  if (!inherits(object, "markovchain")) {
    stop("object must be a 'markovchain' object", call. = FALSE)
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot.markovchain()", call. = FALSE)
  }

  type <- match.arg(type)
  if (type == "eigenvalues") {
    return(.autoplotEigenvalues(object))
  }
  if (type == "flow") {
    return(.autoplotFlow(object, steps = steps, initial = initial))
  }
  if (type == "comparison") {
    return(.autoplotComparison(object, other = other, what = match.arg(what),
                               digits = digits,
                               show_probabilities = show_probabilities))
  }

  if (length(threshold) != 1L || !is.numeric(threshold) ||
      is.na(threshold) || threshold < 0 || threshold > 1) {
    stop("threshold must be a single number between 0 and 1", call. = FALSE)
  }

  if (length(digits) != 1L || !is.numeric(digits) || is.na(digits) ||
      digits < 0) {
    stop("digits must be a non-negative number", call. = FALSE)
  }

  mat <- object@transitionMatrix
  if (!object@byrow) {
    mat <- t(mat)
  }

  states <- object@states
  n <- length(states)
  node_radius <- 0.13

  theta <- seq(0, 2 * pi, length.out = n + 1L)[-(n + 1L)] + pi / 2
  nodes <- data.frame(
    state = states,
    x = cos(theta),
    y = sin(theta),
    stringsAsFactors = FALSE
  )

  classes <- communicatingClasses(object)
  class_id <- integer(n)
  for (i in seq_along(classes)) {
    class_id[match(classes[[i]], states)] <- i
  }
  nodes$class <- factor(class_id)

  edge_index <- which(mat > threshold, arr.ind = TRUE)
  edges <- data.frame()
  loops <- data.frame()

  if (nrow(edge_index) > 0L) {
    edge_rows <- vector("list", nrow(edge_index))
    loop_rows <- list()

    for (k in seq_len(nrow(edge_index))) {
      i <- edge_index[k, 1]
      j <- edge_index[k, 2]
      probability <- mat[i, j]

      if (i == j) {
        a <- seq(-pi / 2, 3 * pi / 2, length.out = 40)
        loop_rows[[length(loop_rows) + 1L]] <- data.frame(
          x = nodes$x[i] + 0.10 * cos(a),
          y = nodes$y[i] + 0.10 + 0.10 * sin(a),
          group = paste0("loop_", i),
          probability = probability,
          stringsAsFactors = FALSE
        )
      } else {
        dx <- nodes$x[j] - nodes$x[i]
        dy <- nodes$y[j] - nodes$y[i]
        distance <- sqrt(dx^2 + dy^2)
        x_start <- nodes$x[i] + node_radius * dx / distance
        y_start <- nodes$y[i] + node_radius * dy / distance
        x_end <- nodes$x[j] - node_radius * dx / distance
        y_end <- nodes$y[j] - node_radius * dy / distance

        edge_rows[[k]] <- data.frame(
          x = x_start,
          y = y_start,
          xend = x_end,
          yend = y_end,
          probability = probability,
          label_x = (x_start + x_end) / 2,
          label_y = (y_start + y_end) / 2,
          stringsAsFactors = FALSE
        )
      }
    }

    edge_rows <- edge_rows[!vapply(edge_rows, is.null, logical(1))]
    if (length(edge_rows) > 0L) {
      edges <- do.call(rbind, edge_rows)
    }
    if (length(loop_rows) > 0L) {
      loops <- do.call(rbind, loop_rows)
    }
  }

  p <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::coord_equal(xlim = c(-1.25, 1.25), ylim = c(-1.35, 1.25),
                         expand = FALSE)

  if (nrow(edges) > 0L) {
    p <- p + ggplot2::geom_curve(
      data = edges,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                   linewidth = probability),
      curvature = 0.15,
      lineend = "round",
      arrow = grid::arrow(length = grid::unit(0.14, "cm"), type = "closed")
    ) +
      ggplot2::scale_linewidth(
        range = c(edge_width * 0.6, edge_width * 1.6),
        guide = "none"
      )

    if (isTRUE(show_probabilities)) {
      p <- p + ggplot2::geom_label(
        data = edges,
        ggplot2::aes(x = label_x, y = label_y,
                     label = formatC(probability, format = "f", digits = digits)),
        size = 3,
        linewidth = 0,
        fill = "white"
      )
    }
  }

  if (nrow(loops) > 0L) {
    p <- p + ggplot2::geom_path(
      data = loops,
      ggplot2::aes(x = x, y = y, group = group, linewidth = probability),
      lineend = "round",
      arrow = grid::arrow(length = grid::unit(0.14, "cm"), type = "closed")
    ) +
      ggplot2::scale_linewidth(
        range = c(edge_width * 0.6, edge_width * 1.6),
        guide = "none"
      )

    if (isTRUE(show_probabilities)) {
      loop_labels <- unique(loops[c("group", "probability")])
      loop_labels$x <- nodes$x[match(sub("loop_", "", loop_labels$group),
                                     seq_len(n))]
      loop_labels$y <- nodes$y[match(sub("loop_", "", loop_labels$group),
                                     seq_len(n))] + 0.25
      p <- p + ggplot2::geom_label(
        data = loop_labels,
        ggplot2::aes(x = x, y = y,
                     label = formatC(probability, format = "f", digits = digits)),
        size = 3,
        linewidth = 0,
        fill = "white"
      )
    }
  }

  p +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = x, y = y, fill = class),
      shape = 21,
      size = node_size,
      stroke = 1
    ) +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = x, y = y, label = state),
      size = 3.5,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_discrete(name = "Communicating class") +
    ggplot2::labs(title = object@name) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}

# Eigenvalues of the transition matrix in the complex plane.
.autoplotEigenvalues <- function(object) {
  P <- .rowStochasticMatrix(object)
  values <- eigen(P, only.values = TRUE)$values
  ev <- data.frame(re = Re(values), im = Im(values))
  theta <- seq(0, 2 * pi, length.out = 361L)
  circle <- data.frame(x = cos(theta), y = sin(theta))

  ggplot2::ggplot() +
    ggplot2::geom_path(data = circle, ggplot2::aes(x = x, y = y),
                       linetype = "dashed", colour = "grey50") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
    ggplot2::geom_point(data = ev, ggplot2::aes(x = re, y = im),
                        size = 3, shape = 21, fill = "steelblue") +
    ggplot2::coord_equal(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) +
    ggplot2::labs(title = object@name,
                  subtitle = "Eigenvalues of the transition matrix",
                  x = "Real part", y = "Imaginary part") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5))
}

# Evolution of the distribution over the states.
.autoplotFlow <- function(object, steps, initial) {
  traj <- redistribute(object, steps = steps, initial = initial)
  long <- data.frame(
    step = rep(as.numeric(rownames(traj)), times = ncol(traj)),
    state = factor(rep(colnames(traj), each = nrow(traj)),
                   levels = colnames(traj)),
    probability = as.vector(traj)
  )

  ggplot2::ggplot(long, ggplot2::aes(x = step, y = probability,
                                     colour = state, group = state)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = object@name,
                  subtitle = "Evolution of the distribution",
                  x = "Step", y = "Probability", colour = "State") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5))
}

# Comparison of two or more chains defined on the same states.
.autoplotComparison <- function(object, other, what, digits,
                                show_probabilities) {
  if (is.null(other)) {
    stop("'other' must be supplied when type = \"comparison\"", call. = FALSE)
  }
  if (inherits(other, "markovchain")) {
    other <- list(other)
  }
  if (!is.list(other) || length(other) == 0L ||
      !all(vapply(other, inherits, logical(1), what = "markovchain"))) {
    stop("'other' must be a markovchain object or a list of markovchain objects",
         call. = FALSE)
  }

  chains <- c(list(object), unname(other))
  labels <- c(object@name,
              if (is.null(names(other))) vapply(other, function(x) x@name, "")
              else names(other))
  labels[is.na(labels) | !nzchar(labels)] <- paste("Chain",
    which(is.na(labels) | !nzchar(labels)))
  labels <- make.unique(labels, sep = " ")

  states <- object@states
  mats <- lapply(chains, function(ch) {
    if (!setequal(ch@states, states) || anyDuplicated(ch@states)) {
      stop("all chains must be defined on the same set of states",
           call. = FALSE)
    }
    .rowStochasticMatrix(ch)[states, states, drop = FALSE]
  })

  if (what == "stationary") {
    pis <- lapply(seq_along(chains), function(i) {
      if (!is.irreducible(chains[[i]])) {
        stop("what = \"stationary\" requires irreducible chains; '",
             labels[i], "' is not irreducible", call. = FALSE)
      }
      dist <- steadyStates(chains[[i]])[1L, ]
      as.numeric(dist[states])
    })
    long <- data.frame(
      chain = factor(rep(labels, each = length(states)), levels = labels),
      state = factor(rep(states, times = length(chains)), levels = states),
      probability = unlist(pis)
    )
    return(
      ggplot2::ggplot(long, ggplot2::aes(x = state, y = probability,
                                         fill = chain)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                          width = 0.7) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::labs(title = "Stationary distributions",
                      x = "State", y = "Probability", fill = "Chain") +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    )
  }

  n <- length(states)
  long <- do.call(rbind, lapply(seq_along(mats), function(i) {
    data.frame(
      chain = factor(labels[i], levels = labels),
      from = factor(rep(states, times = n), levels = rev(states)),
      to = factor(rep(states, each = n), levels = states),
      probability = as.vector(mats[[i]])
    )
  }))
  long$label <- formatC(long$probability, format = "f", digits = digits)

  p <- ggplot2::ggplot(long, ggplot2::aes(x = to, y = from,
                                          fill = probability)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue",
                                 limits = c(0, 1), name = "Probability") +
    ggplot2::facet_wrap(~chain) +
    ggplot2::labs(title = "Transition matrices", x = "To", y = "From") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   panel.grid = ggplot2::element_blank())
  if (isTRUE(show_probabilities)) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = label), size = 3)
  }
  p
}
