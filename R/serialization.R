# Overview: a chain is saved/loaded in two layers.
#   1. toDictionary()/fromDictionary(): chain <-> nested R list (name, states, byrow, matrix as list of lists);
#   2. toFile()/fromFile(): the dictionary <-> a file in JSON, YAML, CSV or XML.
# The JSON/YAML/XML packages are optional (Suggests), so they are only used when needed.
# Internal helper: build the states x states nested named-list representation
# shared by toDictionary() and every toFile() format, from a row-stochastic
# matrix with dimnames already set to the state names.
.matrixToNestedList <- function(P, stateNames) {
  # For every departure state build a named list of probabilities towards all states.
  stats::setNames(lapply(stateNames, function(from) {
    stats::setNames(as.list(as.numeric(P[from, ])), stateNames)
  }), stateNames)
}

# Internal helper: the inverse of .matrixToNestedList(), also accepting a
# plain n x n matrix (with or without dimnames) for convenience when a
# dictionary is built by hand rather than round-tripped through toDictionary().
# The reverse operation: nested list (or plain matrix) -> numeric matrix with state names.
.nestedListToMatrix <- function(tm, stateNames) {
  n <- length(stateNames)
  P <- matrix(NA_real_, n, n, dimnames = list(stateNames, stateNames))

  # Case A: a plain matrix. It must be n x n; if it has names they must match the states.
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
  # Case B: a nested list tm[[from]][[to]]; check that every state appears as a key, then fill row by row.
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
# Always export row-stochastic: transpose if the chain is stored by column.
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
# Rebuilds a chain from a dictionary, validating the content before calling new().
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
# Choose the file format from the extension (.yml is an alias of .yaml), or validate the explicit one.
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
# Write to file: first make the dictionary, then write it in the requested format.
setMethod("toFile", "markovchain", function(object, file, format = NULL) {
  if (length(file) != 1L || !is.character(file) || !nzchar(file)) {
    stop("file must be a single non-empty file path.")
  }
  format <- .resolveSerializationFormat(file, format)
  d <- toDictionary(object)

  if (format == "json") {
    # JSON and YAML need optional packages: stop with a clear message if they are missing.
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
    # CSV: only the matrix, with the state names as header row and first column.
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
# Read from file: the format decides how to read; JSON/YAML/XML go through the dictionary or the
# XML reader, CSV is read directly.
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
    d <- yaml::read_yaml(file, eval.expr = FALSE)
    return(fromDictionary(d))
  }

  if (format == "xml") {
    return(.readChainXml(file))
  }

  # CSV reading: the header row and the first column must list the same states in the same order.
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


# Escapes the five special XML characters (& < > " ') so names can be written safely as attributes.
# Internal helper: escape a string for use inside a double-quoted XML
# attribute.
.xmlEscape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  gsub("'", "&apos;", x, fixed = TRUE)
}

# Writes the PyDTMC XML format: one <Item state_from= state_to= probability=/> for each pair.
# Internal helper: write the dictionary of a chain in PyDTMC's XML format,
# plus the name as an attribute of the root element. Built as text, so that
# writing needs no XML package.
.writeChainXml <- function(d, file) {
  states <- d$states
  from <- rep(states, each = length(states))
  to <- rep(states, times = length(states))
  probs <- unlist(lapply(states, function(s)
    unlist(d$transitionMatrix[[s]][states], use.names = FALSE)))
  # sprintf('%.17g'): 17 significant digits are enough to recover any double exactly.
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

# Reads decimal strings with EXACT rounding. Why? R's as.numeric() can be one unit in the last digit
# away from the closest double on some platforms (e.g. Apple silicon), which would break the exact
# round trip toFile() -> fromFile(). The code below implements 'error-free' arithmetic where a number
# is stored as a sum of two doubles (hi + lo), so no precision is lost before the final rounding.
# Junior note: you do not need to modify this helper; it is covered by tests.
# Internal helper: convert decimal strings to doubles, correctly rounded.
#
# as.numeric() uses R's own parser, which accumulates the digits in a
# `long double`. On platforms where `long double` is just a `double` (notably
# arm64 macOS) a string with 16-17 significant digits can then come out one
# unit in the last place away from the nearest double, so the 17 digits that
# toFile() writes would not bring a probability back exactly. Here the integer
# made of the digits is held as an unevaluated sum of two doubles (exact up to
# 17 digits), divided or multiplied by powers of ten that are exact doubles
# (10^22 and 10^r, r < 22) one after the other, and rounded once: the result
# does not depend on the precision of `long double`. Strings outside this range
# (more than 17 significant digits, a decimal exponent beyond 280 in absolute
# value, or not a plain decimal number) fall back to as.numeric().
.asDoubleExact <- function(s) {
  # v = the fallback result, R's own parser; it is kept for every string the code below
  # does not handle (and for NA)
  v <- suppressWarnings(as.numeric(s))
  s <- trimws(as.character(s))
  # only plain decimal numbers (optional sign, digits with an optional point, optional exponent)
  plain <- !is.na(s) & grepl("^[+-]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][+-]?[0-9]+)?$", s)
  if (!any(plain)) {
    return(v)
  }
  idx <- which(plain)
  x <- s[idx]
  # Split each number as sign, decimal exponent `ex`, and the digits before and after the point.
  neg <- startsWith(x, "-")
  x <- sub("^[+-]", "", x)
  ex <- integer(length(x))
  hasExp <- grepl("[eE]", x)
  ex[hasExp] <- suppressWarnings(as.integer(sub("^.*[eE]", "", x[hasExp])))
  x <- sub("[eE].*$", "", x)
  int <- sub("\\..*$", "", x)
  frac <- ifelse(grepl("\\.", x), sub("^[^.]*\\.", "", x), "")
  # The value is then the integer `digits` times 10^e10 (leading zeros removed).
  digits <- sub("^0+", "", paste0(int, frac))
  nd <- nchar(digits)
  e10 <- ex - nchar(frac)
  # handled here only if the integer has at most 17 digits and the power of ten is moderate
  ok <- !is.na(e10) & nd >= 1L & nd <= 17L & abs(e10) <= 280L
  if (!any(ok)) {
    return(v)
  }
  k <- which(ok)
  digits <- digits[k]
  nd <- nd[k]
  e10 <- e10[k]
  # The integer made of all the digits of the number is split in two parts that fit exactly in a double.
  # the integer as hi * 1e8 + lo, both parts exact doubles (hi < 1e9, lo < 1e8)
  hi <- ifelse(nd > 8L, as.numeric(substr(digits, 1L, nd - 8L)), 0)
  lo <- as.numeric(substr(digits, pmax(nd - 7L, 1L), nd))
  # twoProd/twoSum: classic algorithms (Dekker/Veltkamp, Knuth) returning a product/sum together with
  # its exact rounding error.
  # error-free transformations (Veltkamp-Dekker), vectorised
  twoProd <- function(a, b) {
    p <- a * b
    ca <- 134217729 * a
    ah <- ca - (ca - a)
    al <- a - ah
    cb <- 134217729 * b
    bh <- cb - (cb - b)
    bl <- b - bh
    list(p = p, e = ((ah * bh - p) + ah * bl + al * bh) + al * bl)
  }
  twoSum <- function(a, b) {
    s <- a + b
    bb <- s - a
    list(s = s, e = (a - (s - bb)) + (b - bb))
  }
  # (nhi, nlo) = hi * 1e8 + lo as an unevaluated sum of two doubles, exact for up to 17 digits
  pr <- twoProd(hi, 1e8)
  sm <- twoSum(pr$p, lo)
  lo2 <- sm$e + pr$e
  nhi <- sm$s + lo2
  nlo <- lo2 - (nhi - sm$s)
  # (hi, lo) divided or multiplied by an exact double d, keeping both parts
  divide <- function(hi, lo, d) {
    q1 <- hi / d
    pp <- twoProd(q1, d)
    q2 <- (((hi - pp$p) - pp$e) + lo) / d
    s <- q1 + q2
    list(hi = s, lo = q2 - (s - q1))
  }
  multiply <- function(hi, lo, d) {
    pp <- twoProd(hi, d)
    l <- pp$e + lo * d
    s <- pp$p + l
    list(hi = s, lo = l - (s - pp$p))
  }
  # The power of ten is applied in steps of 1e22 (the largest exact power of ten in a double), then a remainder.
  # |e10| = 22 * nBig + r: nBig steps by 1e22, then one by 10^r (all exact)
  nBig <- abs(e10) %/% 22L
  pw <- cumprod(c(1, rep(10, 21L)))[abs(e10) %% 22L + 1L]
  down <- e10 < 0L                      # negative exponent: divide; otherwise multiply
  cur <- list(hi = nhi, lo = nlo)
  # step(): applies one division or multiplication by d to the elements selected by `sel`
  step <- function(cur, sel, d) {
    if (!any(sel)) {
      return(cur)
    }
    out <- list(hi = cur$hi, lo = cur$lo)
    for (dir in c(TRUE, FALSE)) {
      m <- sel & (down == dir)
      if (any(m)) {
        res <- (if (dir) divide else multiply)(cur$hi[m], cur$lo[m], if (length(d) == 1L) d else d[m])
        out$hi[m] <- res$hi
        out$lo[m] <- res$lo
      }
    }
    out
  }
  for (i in seq_len(max(nBig))) {
    cur <- step(cur, nBig >= i, 1e22)
  }
  cur <- step(cur, pw != 1, pw)
  # the single final rounding hi + lo; put the results back in their positions, keeping the
  # fallback value of as.numeric() where the exact computation did not give a finite number
  pos <- idx[k]
  res <- cur$hi + cur$lo
  res <- ifelse(neg[k], -res, res)
  good <- is.finite(res)
  v[pos[good]] <- res[good]
  v
}

# Reads PyDTMC's XML: checks the structure, collects (from, to, probability) triples, checks that all
# n*n pairs are present exactly once, and fills the matrix.
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
  # each <Item state_from="a" state_to="b" probability="0.3"/> is one cell of the matrix
  required <- c("probability", "state_from", "state_to")
  attrs <- xml2::xml_attrs(items)
  if (any(vapply(attrs, function(a) !identical(sort(names(a)), required), logical(1)))) {
    stop("Every 'Item' element must have exactly the attributes state_from, state_to and probability.")
  }
  from <- trimws(vapply(attrs, `[[`, character(1), "state_from"))
  to <- trimws(vapply(attrs, `[[`, character(1), "state_to"))
  probs <- .asDoubleExact(vapply(attrs, `[[`, character(1), "probability"))
  if (any(!nzchar(from)) || any(!nzchar(to)) || anyNA(probs)) {
    stop("The XML file contains empty state names or probabilities that are not numbers.")
  }
  # The file lists every pair (from, to), self transitions included, so the states are
  # those that appear as a self transition; there must be exactly n * n distinct items.
  stateNames <- unique(from[from == to])
  n <- length(stateNames)
  if (n == 0L || !setequal(unique(c(from, to)), stateNames) ||
      length(items) != n * n || anyDuplicated(paste(from, to, sep = "\r"))) {
    stop("The XML file must contain exactly one 'Item' for every pair of states, self transitions included.")
  }
  P <- matrix(0, n, n, dimnames = list(stateNames, stateNames))
  P[cbind(match(from, stateNames), match(to, stateNames))] <- probs   # [from, to]
  nm <- xml2::xml_attr(doc, "name")
  if (is.na(nm)) {
    nm <- "Unnamed Markov chain"
  }
  new("markovchain", states = stateNames, byrow = TRUE, transitionMatrix = P, name = nm)
}
