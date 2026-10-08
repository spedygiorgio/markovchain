# Internal helper: build the states x states nested named-list representation
# shared by toDictionary() and every toFile() format, from a row-stochastic
# matrix with dimnames already set to the state names.
.matrixToNestedList <- function(P, stateNames) {
  stats::setNames(lapply(stateNames, function(from) {
    stats::setNames(as.list(as.numeric(P[from, ])), stateNames)
  }), stateNames)
}

# Internal helper: the inverse of .matrixToNestedList(), also accepting a
# plain n x n matrix (with or without dimnames) for convenience when a
# dictionary is built by hand rather than round-tripped through toDictionary().
.nestedListToMatrix <- function(tm, stateNames) {
  n <- length(stateNames)
  P <- matrix(NA_real_, n, n, dimnames = list(stateNames, stateNames))

  if (is.matrix(tm)) {
    if (!identical(dim(tm), c(n, n))) {
      stop("transitionMatrix must be an n x n matrix, with n = length(states).")
    }
    rn <- rownames(tm)
    cn <- colnames(tm)
    if (!is.null(rn) || !is.null(cn)) {
      if (is.null(rn) || is.null(cn) || !setequal(rn, stateNames) || !setequal(cn, stateNames)) {
        stop("transitionMatrix's row/column names must match states exactly.")
      }
      P <- tm[stateNames, stateNames, drop = FALSE]
    } else {
      P[] <- as.numeric(tm)
    }
  } else if (is.list(tm)) {
    if (is.null(names(tm)) || !setequal(names(tm), stateNames)) {
      stop("transitionMatrix must be named after every state (one entry per source state).")
    }
    for (from in stateNames) {
      thisRow <- tm[[from]]
      if (is.null(thisRow) || is.null(names(thisRow)) || !setequal(names(thisRow), stateNames)) {
        stop(paste0("transitionMatrix[[\"", from, "\"]] must be a named list/vector with one entry per state."))
      }
      P[from, stateNames] <- as.numeric(unlist(thisRow[stateNames], use.names = FALSE))
    }
  } else {
    stop("transitionMatrix must be a matrix or a named list of named lists/vectors (as produced by toDictionary()).")
  }

  storage.mode(P) <- "double"
  P
}

#' Represent a Markov chain as a plain R list
#'
#' Converts a \code{markovchain} object to a plain, self-describing R list:
#' the same information \code{\link{toFile}} writes to disk, kept in memory.
#' \code{fromDictionary} reverses the conversion.
#'
#' @param object A \code{markovchain} object.
#'
#' @return A named list with four elements:
#'   \describe{
#'     \item{\code{name}}{The chain's \code{name}, as a single string
#'       (possibly empty).}
#'     \item{\code{states}}{A character vector of state names, in order.}
#'     \item{\code{byrow}}{Always \code{TRUE}: the list always stores the
#'       chain row-stochastically, regardless of \code{object}'s own
#'       storage convention, so that the representation is unambiguous
#'       without also having to interpret this flag.}
#'     \item{\code{transitionMatrix}}{A named list of named lists:
#'       \code{transitionMatrix[[i]][[j]]} is the probability of moving
#'       from state \code{i} to state \code{j}. This is deliberately not a
#'       plain matrix, so that the structure serializes to JSON or YAML
#'       (via \code{\link{toFile}}) as a self-describing object keyed by
#'       state name, rather than a bare array whose meaning depends on
#'       remembering a row/column order.}
#'   }
#'
#' @details
#' Unlike PyDTMC's own \code{to_dictionary()}/\code{from_dictionary()},
#' which represent a chain as a flat mapping from every
#' \code{(from_state, to_state)} pair to its probability -- \eqn{n^2}
#' entries with no state grouping -- this nests the representation by
#' source state, which is both more compact to read and directly
#' round-trips through R's own list-of-lists idiom without any special
#' tuple-key handling.
#'
#' @seealso \code{\link{toFile}}, \code{\link{fromFile}}
#'
#' @examples
#' statesNames <- c("a", "b")
#' mc <- new("markovchain", states = statesNames,
#'           transitionMatrix = matrix(c(0.7, 0.3, 0.4, 0.6), byrow = TRUE,
#'                                      nrow = 2, dimnames = list(statesNames, statesNames)))
#' d <- toDictionary(mc)
#' d$transitionMatrix$a$b # 0.3: probability of moving from "a" to "b"
#'
#' identical(fromDictionary(d)@transitionMatrix, mc@transitionMatrix)
#'
#' @exportMethod toDictionary
setGeneric("toDictionary", function(object) standardGeneric("toDictionary"))

#' @rdname toDictionary
setMethod("toDictionary", "markovchain", function(object) {
  P <- as.matrix(object@transitionMatrix)
  if (!object@byrow) {
    P <- t(P)
  }
  stateNames <- states(object)
  dimnames(P) <- list(stateNames, stateNames)

  list(
    name = object@name,
    states = stateNames,
    byrow = TRUE,
    transitionMatrix = .matrixToNestedList(P, stateNames)
  )
})

#' @rdname toDictionary
#'
#' @param d A list as returned by \code{toDictionary}. \code{d$transitionMatrix}
#'   may also be a plain \code{n x n} matrix (with or without dimnames
#'   matching \code{d$states}) for convenience when building a dictionary by
#'   hand rather than from an existing \code{markovchain} object. When it is
#'   a plain matrix, \code{d$byrow} is honored exactly as
#'   \code{\link{new}("markovchain", ...)} honors its own \code{byrow}
#'   argument: the matrix is stored as given, with \code{d$byrow} only
#'   documenting whether it is row- or column-stochastic (so a
#'   column-stochastic matrix round-trips by setting \code{d$byrow = FALSE},
#'   with no transposition performed here). When \code{d$transitionMatrix}
#'   is a nested list (as \code{toDictionary} produces), the nesting itself
#'   is always keyed \code{[[from]][[to]]} -- i.e. row-stochastic -- so it is
#'   always reconstructed with \code{byrow = TRUE}, regardless of
#'   \code{d$byrow}.
#' @return \code{fromDictionary} returns a \code{markovchain} object.
#' @export
fromDictionary <- function(d) {
  if (!is.list(d)) {
    stop("d must be a list, as returned by toDictionary().")
  }
  if (is.null(d$states) || !is.character(d$states) || length(d$states) < 1L ||
      anyNA(d$states) || anyDuplicated(d$states)) {
    stop("d$states must be a character vector of unique, non-missing state names.")
  }
  if (is.null(d$transitionMatrix)) {
    stop("d$transitionMatrix is missing.")
  }

  stateNames <- d$states
  P <- .nestedListToMatrix(d$transitionMatrix, stateNames)

  # A plain matrix carries its own row/column-stochastic convention, exactly
  # like new("markovchain", ...): honor d$byrow without transposing. A
  # nested list is always keyed [[from]][[to]] (row-stochastic by
  # construction), so it is always reconstructed with byrow = TRUE.
  byrow <- if (is.matrix(d$transitionMatrix)) {
    if (is.null(d$byrow)) TRUE else isTRUE(d$byrow)
  } else {
    TRUE
  }

  nm <- if (is.null(d$name) || length(d$name) != 1L || is.na(d$name)) {
    "Unnamed Markov chain" # same default new("markovchain", ...) itself uses
  } else {
    as.character(d$name)
  }

  new("markovchain", states = stateNames, byrow = byrow, transitionMatrix = P, name = nm)
}

# Internal helper: infer a supported serialization format from a file's
# extension, or validate an explicitly given one.
.resolveSerializationFormat <- function(file, format) {
  supported <- c("json", "yaml", "csv", "xml")
  if (is.null(format)) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "yml") {
      ext <- "yaml"
    }
    if (!(ext %in% supported)) {
      stop(paste0(
        "Unable to infer a format from the file extension of \"", file,
        "\". Pass format explicitly: one of \"json\", \"yaml\", \"csv\" or \"xml\"."
      ))
    }
    return(ext)
  }
  format <- match.arg(format, supported)
  format
}

#' Write or read a Markov chain to or from a file
#'
#' Writes a \code{markovchain} object to a JSON, YAML, CSV or XML file, or reads
#' one back, using the same representation as \code{\link{toDictionary}}.
#'
#' @param object A \code{markovchain} object (for \code{toFile}).
#' @param file A single file path to write to or read from. If \code{format}
#'   is not supplied, it is inferred from the file extension
#'   (\code{.json}, \code{.yaml}/\code{.yml}, \code{.csv} or \code{.xml}).
#' @param format One of \code{"json"}, \code{"yaml"}, \code{"csv"} or
#'   \code{"xml"}. The
#'   default, \code{NULL}, infers the format from \code{file}'s extension.
#'
#' @return \code{toFile} returns \code{file}, invisibly. \code{fromFile}
#'   returns a \code{markovchain} object.
#'
#' @details
#' The JSON and YAML formats store exactly what \code{\link{toDictionary}}
#' returns (name, state names, and the transition probabilities nested by
#' source state), and round-trip a chain exactly, including full numeric
#' precision (\code{toFile} writes both JSON and YAML with 17 significant
#' digits, enough to recover every double exactly).
#'
#' The CSV format only stores the transition matrix itself, as a table of
#' probabilities with the state names as both the header row and the first
#' column -- there is no natural place in a CSV file for the chain's
#' \code{name}, so it is not preserved by \code{toFile(..., format = "csv")}
#' and \code{fromFile} always returns an unnamed chain for a \code{.csv}
#' file. This is the same limitation PyDTMC's own CSV format has.
#'
#' The XML format is the one of PyDTMC, so files can be exchanged with it in
#' both directions: a root element \code{MarkovChain} with one \code{Item}
#' element per transition, whose attributes are \code{state_from},
#' \code{state_to} and \code{probability}. All \eqn{n^2} transitions are
#' written, zeros included, and probabilities use 17 significant digits, so
#' the round trip is exact. The \code{name} of the chain is stored as an
#' attribute of the root element, which PyDTMC ignores when reading. When
#' reading, the states are taken in the order in which their self
#' transitions (\code{state_from} equal to \code{state_to}) appear, as
#' PyDTMC does, and the name is restored if the attribute is present.
#' Writing XML uses only base R; reading it requires the \pkg{xml2}
#' package.
#'
#' Writing JSON requires the \pkg{jsonlite} package, and writing YAML
#' requires the \pkg{yaml} package; both are only in \code{Suggests}, and an
#' informative error is raised if the relevant package is not installed.
#' Reading has the same requirements for the format being read. CSV uses
#' only base R and has no extra package dependency.
#'
#' @seealso \code{\link{toDictionary}}, \code{\link{fromDictionary}}
#'
#' @examples
#' \dontrun{
#' statesNames <- c("a", "b")
#' mc <- new("markovchain", states = statesNames,
#'           transitionMatrix = matrix(c(0.7, 0.3, 0.4, 0.6), byrow = TRUE,
#'                                      nrow = 2, dimnames = list(statesNames, statesNames)))
#' toFile(mc, "chain.json")
#' identical(fromFile("chain.json")@transitionMatrix, mc@transitionMatrix)
#' }
#'
#' @exportMethod toFile
setGeneric("toFile", function(object, file, format = NULL) standardGeneric("toFile"))

#' @rdname toFile
setMethod("toFile", "markovchain", function(object, file, format = NULL) {
  if (length(file) != 1L || !is.character(file) || !nzchar(file)) {
    stop("file must be a single non-empty file path.")
  }
  format <- .resolveSerializationFormat(file, format)
  d <- toDictionary(object)

  if (format == "json") {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Writing JSON requires the 'jsonlite' package. Install it with install.packages(\"jsonlite\").")
    }
    jsonlite::write_json(d, file, auto_unbox = TRUE, pretty = TRUE, digits = 17)
  } else if (format == "yaml") {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Writing YAML requires the 'yaml' package. Install it with install.packages(\"yaml\").")
    }
    yaml::write_yaml(d, file, precision = 17)
  } else if (format == "csv") {
    P <- matrix(unlist(lapply(d$states, function(from) unlist(d$transitionMatrix[[from]][d$states], use.names = FALSE))),
                nrow = length(d$states), byrow = TRUE,
                dimnames = list(d$states, d$states))
    utils::write.csv(P, file, row.names = TRUE)
  } else {
    .writeChainXml(d, file)
  }

  invisible(file)
})

#' @rdname toFile
#' @export
fromFile <- function(file, format = NULL) {
  if (length(file) != 1L || !is.character(file) || !nzchar(file)) {
    stop("file must be a single non-empty file path.")
  }
  if (!file.exists(file)) {
    stop(paste0("File not found: ", file))
  }
  format <- .resolveSerializationFormat(file, format)

  if (format == "json") {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Reading JSON requires the 'jsonlite' package. Install it with install.packages(\"jsonlite\").")
    }
    d <- jsonlite::read_json(file, simplifyVector = TRUE)
    return(fromDictionary(d))
  }
  if (format == "yaml") {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Reading YAML requires the 'yaml' package. Install it with install.packages(\"yaml\").")
    }
    d <- yaml::read_yaml(file)
    return(fromDictionary(d))
  }

  if (format == "xml") {
    return(.readChainXml(file))
  }

  # CSV: the header row and the first column both give the state names; no
  # "name" field is stored (see Details in ?toFile), so the result is unnamed.
  raw <- utils::read.csv(file, row.names = 1, check.names = FALSE)
  stateNames <- rownames(raw)
  if (!identical(colnames(raw), stateNames)) {
    stop("The CSV file's header row and first column must list the same state names, in the same order.")
  }
  P <- as.matrix(raw)
  dimnames(P) <- list(stateNames, stateNames)
  new("markovchain", states = stateNames, byrow = TRUE, transitionMatrix = P)
}


# Internal helper: escape a string for use inside a double-quoted XML
# attribute.
.xmlEscape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  gsub("'", "&apos;", x, fixed = TRUE)
}

# Internal helper: write the dictionary of a chain in PyDTMC's XML format,
# plus the name as an attribute of the root element. Built as text, so that
# writing needs no XML package.
.writeChainXml <- function(d, file) {
  states <- d$states
  from <- rep(states, each = length(states))
  to <- rep(states, times = length(states))
  probs <- unlist(lapply(states, function(s)
    unlist(d$transitionMatrix[[s]][states], use.names = FALSE)))
  items <- paste0("\t<Item state_from=\"", .xmlEscape(from),
                  "\" state_to=\"", .xmlEscape(to),
                  "\" probability=\"", sprintf("%.17g", probs), "\"/>")
  lines <- c("<?xml version='1.0' encoding='utf-8' standalone='yes' ?>",
             paste0("<MarkovChain name=\"", .xmlEscape(d$name), "\">"),
             items,
             "</MarkovChain>")
  con <- file(file, open = "w", encoding = "UTF-8")
  on.exit(close(con))
  writeLines(lines, con)
}

# Internal helper: read a chain from PyDTMC's XML format.
.readChainXml <- function(file) {
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("Reading XML requires the 'xml2' package. Install it with install.packages(\"xml2\").")
  }
  doc <- tryCatch(xml2::read_xml(file),
                  error = function(e) stop("The XML file could not be parsed: ", conditionMessage(e)))
  if (xml2::xml_name(doc) != "MarkovChain") {
    stop("The root element of the XML file must be 'MarkovChain'.")
  }
  items <- xml2::xml_children(doc)
  if (length(items) == 0L || any(xml2::xml_name(items) != "Item")) {
    stop("The XML file must contain only 'Item' elements.")
  }
  required <- c("probability", "state_from", "state_to")
  attrs <- xml2::xml_attrs(items)
  if (any(vapply(attrs, function(a) !identical(sort(names(a)), required), logical(1)))) {
    stop("Every 'Item' element must have exactly the attributes state_from, state_to and probability.")
  }
  from <- trimws(vapply(attrs, `[[`, character(1), "state_from"))
  to <- trimws(vapply(attrs, `[[`, character(1), "state_to"))
  probs <- suppressWarnings(as.numeric(vapply(attrs, `[[`, character(1), "probability")))
  if (any(!nzchar(from)) || any(!nzchar(to)) || anyNA(probs)) {
    stop("The XML file contains empty state names or probabilities that are not numbers.")
  }
  stateNames <- unique(from[from == to])
  n <- length(stateNames)
  if (n == 0L || !setequal(unique(c(from, to)), stateNames) ||
      length(items) != n * n || anyDuplicated(paste(from, to, sep = "\r"))) {
    stop("The XML file must contain exactly one 'Item' for every pair of states, self transitions included.")
  }
  P <- matrix(0, n, n, dimnames = list(stateNames, stateNames))
  P[cbind(match(from, stateNames), match(to, stateNames))] <- probs
  nm <- xml2::xml_attr(doc, "name")
  if (is.na(nm)) {
    nm <- "Unnamed Markov chain"
  }
  new("markovchain", states = stateNames, byrow = TRUE, transitionMatrix = P, name = nm)
}
