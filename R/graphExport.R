# Export of the transition graph of a Markov chain as Graphviz DOT or
# Mermaid text, for rendering outside R (Graphviz, documentation sites,
# Markdown files that support Mermaid diagrams).

# Internal helper: the edges of the transition graph, one row per transition
# with probability above `minProbability`, in row-major state order.
.chainEdges <- function(object, minProbability) {
  if (length(minProbability) != 1L || !is.numeric(minProbability) ||
      !is.finite(minProbability) || minProbability < 0 ||
      minProbability >= 1) {
    stop("minProbability must be a single number in [0, 1).")
  }
  P <- .rowStochasticMatrix(object)
  idx <- which(t(P) > minProbability, arr.ind = TRUE)
  # t(P) makes which() run along the rows of P: column 1 is the target
  # state, column 2 the source state
  idx <- idx[order(idx[, 2L], idx[, 1L]), , drop = FALSE]
  data.frame(from = idx[, 2L], to = idx[, 1L],
             probability = P[cbind(idx[, 2L], idx[, 1L])])
}

# Internal helper: probability labels with `digits` significant digits.
.edgeLabels <- function(p, digits) {
  if (length(digits) != 1L || !is.numeric(digits) || !is.finite(digits) ||
      digits < 1 || digits != round(digits)) {
    stop("digits must be a single positive whole number.")
  }
  format(signif(p, digits), digits = digits, trim = TRUE, scientific = FALSE,
         drop0trailing = TRUE)
}

# Internal helper: return the text, also writing it to `file` when given.
.graphOutput <- function(lines, file) {
  text <- paste(lines, collapse = "\n")
  if (is.null(file)) {
    return(text)
  }
  if (length(file) != 1L || !is.character(file) || !nzchar(file)) {
    stop("file must be NULL or a single non-empty file path.")
  }
  con <- base::file(file, open = "w", encoding = "UTF-8")
  on.exit(close(con))
  writeLines(lines, con)
  invisible(text)
}

#' Export the transition graph as Graphviz DOT or Mermaid text
#'
#' \code{toDot} writes the transition graph of a Markov chain in the DOT
#' language of Graphviz; \code{toMermaid} writes it as a Mermaid flowchart.
#' The states are the nodes and every transition with positive probability
#' is an edge labelled with its probability.
#'
#' @param object A \code{markovchain} object.
#' @param file An optional file path. If given, the text is also written
#'   there (UTF-8) and returned invisibly.
#' @param digits The number of significant digits of the probabilities shown
#'   on the edges.
#' @param minProbability Transitions with probability not above this value
#'   are left out; the default, 0, keeps every possible transition.
#' @param direction The layout direction: \code{"LR"} (left to right, the
#'   default), \code{"TB"} (top to bottom), \code{"RL"} or \code{"BT"}.
#'
#' @details
#' The edges follow the outgoing distribution of each state whatever the
#' storage of the transition matrix (\code{byrow}), so a chain and its
#' column-stochastic copy give the same text. State names are quoted and
#' escaped, so they may contain spaces, quotes or other punctuation. In the
#' Mermaid output the nodes get the identifiers \code{s1}, \code{s2}, ...,
#' with the state names as labels, because Mermaid identifiers cannot
#' contain arbitrary characters.
#'
#' The DOT text can be rendered with Graphviz (for instance
#' \code{dot -Tpng chain.dot -o chain.png}) or with
#' \code{DiagrammeR::grViz()}; the Mermaid text can be pasted into any
#' Markdown renderer that supports Mermaid diagrams, or rendered with
#' \code{DiagrammeR::mermaid()}.
#'
#' @return The text, as a single character string (invisibly when
#'   \code{file} is given).
#'
#' @seealso \code{\link{toFile}}, \code{\link{plot,markovchain,missing-method}}
#'
#' @examples
#' statesNames <- c("a", "b", "c")
#' mc <- new("markovchain", states = statesNames, name = "Example",
#'   transitionMatrix = matrix(c(0.5, 0.5, 0, 0.2, 0.3, 0.5, 0, 0, 1),
#'     byrow = TRUE, nrow = 3, dimnames = list(statesNames, statesNames)))
#' cat(toDot(mc), "\n")
#' cat(toMermaid(mc, direction = "TB"), "\n")
#'
#' @exportMethod toDot
setGeneric("toDot", function(object, file = NULL, digits = 3,
                             minProbability = 0, direction = "LR")
  standardGeneric("toDot"))

#' @rdname toDot
setMethod("toDot", "markovchain",
  function(object, file = NULL, digits = 3, minProbability = 0,
           direction = "LR") {
    direction <- match.arg(direction, c("LR", "TB", "RL", "BT"))
    edges <- .chainEdges(object, minProbability)
    labels <- .edgeLabels(edges$probability, digits)
    quote <- function(x) {
      x <- gsub("\\", "\\\\", x, fixed = TRUE)
      x <- gsub("\"", "\\\"", x, fixed = TRUE)
      x <- gsub("\n", "\\n", x, fixed = TRUE)
      paste0("\"", x, "\"")
    }
    states <- object@states
    lines <- c(
      paste0("digraph ", quote(object@name), " {"),
      paste0("  rankdir=", direction, ";"),
      "  node [shape=circle];",
      paste0("  ", quote(states), ";"),
      if (nrow(edges) > 0L) {
        paste0("  ", quote(states[edges$from]), " -> ", quote(states[edges$to]),
               " [label=", quote(labels), "];")
      },
      "}"
    )
    .graphOutput(lines, file)
  })

#' @rdname toDot
#' @exportMethod toMermaid
setGeneric("toMermaid", function(object, file = NULL, digits = 3,
                                 minProbability = 0, direction = "LR")
  standardGeneric("toMermaid"))

#' @rdname toDot
setMethod("toMermaid", "markovchain",
  function(object, file = NULL, digits = 3, minProbability = 0,
           direction = "LR") {
    direction <- match.arg(direction, c("LR", "TB", "RL", "BT"))
    edges <- .chainEdges(object, minProbability)
    labels <- .edgeLabels(edges$probability, digits)
    # Mermaid labels in double quotes accept any text except the double
    # quote itself, written as the entity #quot;
    label <- function(x) {
      x <- gsub("\"", "#quot;", x, fixed = TRUE)
      gsub("\n", " ", x, fixed = TRUE)
    }
    states <- object@states
    ids <- paste0("s", seq_along(states))
    lines <- c(
      paste0("flowchart ", direction),
      paste0("  ", ids, "((\"", label(states), "\"))"),
      if (nrow(edges) > 0L) {
        paste0("  ", ids[edges$from], " -->|\"", labels, "\"| ", ids[edges$to])
      }
    )
    .graphOutput(lines, file)
  })
