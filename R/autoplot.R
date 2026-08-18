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
#' }
autoplot.markovchain <- function(object,
                                  threshold = 0,
                                  show_probabilities = TRUE,
                                  digits = 2,
                                  node_size = 6,
                                  edge_width = 1,
                                  ...) {
  if (!inherits(object, "markovchain")) {
    stop("object must be a 'markovchain' object", call. = FALSE)
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot.markovchain()", call. = FALSE)
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
