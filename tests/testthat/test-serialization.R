library(testthat)
library(markovchain)

## ---- fixtures ---------------------------------------------------------

statesRow <- c("a", "b", "c")
Prow <- matrix(c(0.5, 0.3, 0.2,
                  0.2, 0.6, 0.2,
                  0.1, 0.1, 0.8), byrow = TRUE, nrow = 3,
                dimnames = list(statesRow, statesRow))
mcRow <- new("markovchain", states = statesRow, transitionMatrix = Prow, name = "RowChain")

# The transpose of a row-stochastic matrix is column-stochastic, and
# represents the very same chain when interpreted with byrow = FALSE.
Pcol <- t(Prow)
mcCol <- new("markovchain", states = statesRow, byrow = FALSE,
             transitionMatrix = Pcol, name = "ColChain")

# A chain whose entries need full double precision to round-trip exactly.
statesPrec <- c("x", "y")
Pprec <- matrix(c(1 / 3, 2 / 3, 0.123456789012345, 0.876543210987655),
                 byrow = TRUE, nrow = 2, dimnames = list(statesPrec, statesPrec))
mcPrec <- new("markovchain", states = statesPrec, transitionMatrix = Pprec, name = "PrecisionChain")

## ---- toDictionary() / fromDictionary(): row-stochastic round-trip -----

test_that("toDictionary()/fromDictionary() round-trip a row-stochastic chain exactly", {
  d <- toDictionary(mcRow)

  expect_type(d, "list")
  expect_equal(d$name, "RowChain")
  expect_equal(d$states, statesRow)
  expect_true(d$byrow)
  expect_equal(d$transitionMatrix$a$b, unname(Prow["a", "b"]))
  expect_equal(d$transitionMatrix$c$a, unname(Prow["c", "a"]))

  mc2 <- fromDictionary(d)
  expect_s4_class(mc2, "markovchain")
  expect_identical(mc2@transitionMatrix, mcRow@transitionMatrix)
  expect_identical(mc2@name, mcRow@name)
  expect_true(mc2@byrow)
})

## ---- toDictionary() normalizes column-stochastic chains ---------------

test_that("toDictionary() always normalizes to a row-stochastic, byrow = TRUE dictionary", {
  d <- toDictionary(mcCol)

  expect_true(d$byrow)
  # mcCol and mcRow represent the same chain (one is the transpose-encoded
  # form of the other), so their normalized dictionaries must agree exactly.
  expect_equal(d$transitionMatrix, toDictionary(mcRow)$transitionMatrix)

  mc2 <- fromDictionary(d)
  expect_true(mc2@byrow)
  expect_equal(unname(mc2@transitionMatrix), unname(mcRow@transitionMatrix))
})

## ---- fromDictionary() honors byrow for a hand-built plain-matrix dict --

test_that("fromDictionary() honors an explicit byrow for a plain-matrix transitionMatrix, unchanged", {
  # Passing a column-stochastic matrix with byrow = FALSE should reproduce
  # mcCol's own raw (unnormalized) representation exactly, mirroring how
  # new(\"markovchain\", ...) itself treats byrow: no transposition.
  d <- list(states = statesRow, transitionMatrix = Pcol, byrow = FALSE)
  mc2 <- fromDictionary(d)

  expect_false(mc2@byrow)
  expect_identical(unname(mc2@transitionMatrix), unname(mcCol@transitionMatrix))
  expect_equal(unname(colSums(mc2@transitionMatrix)), rep(1, 3))

  # And a plain row-stochastic matrix with the default byrow reproduces mcRow.
  d2 <- list(states = statesRow, transitionMatrix = Prow)
  mc3 <- fromDictionary(d2)
  expect_true(mc3@byrow)
  expect_identical(unname(mc3@transitionMatrix), unname(mcRow@transitionMatrix))
})

## ---- fromDictionary() defaults ------------------------------------------

test_that("fromDictionary() defaults an absent name to the package's own default", {
  d <- list(states = statesRow, transitionMatrix = Prow)
  mc2 <- fromDictionary(d)
  expect_equal(mc2@name, new("markovchain")@name)
})

## ---- fromDictionary() validation ---------------------------------------

test_that("fromDictionary() rejects malformed input", {
  expect_error(fromDictionary("not a list"), "must be a list")
  expect_error(fromDictionary(list(transitionMatrix = Prow)), "states")
  expect_error(fromDictionary(list(states = c("a", "a"), transitionMatrix = Prow)),
               "unique")
  expect_error(fromDictionary(list(states = statesRow)), "transitionMatrix")
  expect_error(
    fromDictionary(list(states = statesRow, transitionMatrix = matrix(1, 2, 2))),
    "n x n matrix"
  )
  expect_error(
    fromDictionary(list(states = statesRow,
                         transitionMatrix = list(a = list(a = 0.5, b = 0.5)))),
    "named after every state"
  )
  expect_error(
    fromDictionary(list(states = statesRow,
                         transitionMatrix = list(a = list(a = 1), b = list(a = 1), c = list(a = 1)))),
    "one entry per state"
  )
  expect_error(fromDictionary(list(states = statesRow, transitionMatrix = 42)),
               "matrix or a named list")
})

## ---- toFile()/fromFile(): format inference ------------------------------

test_that("toFile()/fromFile() infer format from the file extension", {
  skip_if_not_installed("jsonlite")
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp))

  toFile(mcRow, tmp)
  mc2 <- fromFile(tmp)
  expect_identical(mc2@transitionMatrix, mcRow@transitionMatrix)
  expect_identical(mc2@name, mcRow@name)
})

test_that("toFile()/fromFile() reject a file with an unsupported/unknown extension", {
  tmp <- tempfile(fileext = ".toml")  # .xml is supported since 1.2
  expect_error(toFile(mcRow, tmp), "Unable to infer a format")
  expect_error(fromFile(tmp), "not found|Unable to infer a format")
})

test_that("an explicit format overrides the file extension", {
  skip_if_not_installed("jsonlite")
  tmp <- tempfile(fileext = ".dat")
  on.exit(unlink(tmp))

  toFile(mcRow, tmp, format = "json")
  mc2 <- fromFile(tmp, format = "json")
  expect_identical(mc2@transitionMatrix, mcRow@transitionMatrix)
})

test_that("toFile()/fromFile() validate the file argument and file existence", {
  expect_error(toFile(mcRow, character(0)), "single non-empty file path")
  expect_error(toFile(mcRow, c("a.json", "b.json")), "single non-empty file path")
  expect_error(fromFile(tempfile(fileext = ".json")), "not found")
})

## ---- toFile()/fromFile(): JSON round-trip -------------------------------

test_that("JSON round-trips a chain exactly, including full numeric precision", {
  skip_if_not_installed("jsonlite")
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp))

  toFile(mcPrec, tmp)
  mc2 <- fromFile(tmp)
  expect_identical(mc2@transitionMatrix, mcPrec@transitionMatrix)
  expect_identical(mc2@name, mcPrec@name)
  expect_true(mc2@byrow)
})

## ---- toFile()/fromFile(): YAML round-trip -------------------------------

test_that("YAML round-trips a chain exactly, including full numeric precision", {
  skip_if_not_installed("yaml")
  tmp <- tempfile(fileext = ".yaml")
  on.exit(unlink(tmp))

  toFile(mcPrec, tmp)
  mc2 <- fromFile(tmp)
  expect_identical(mc2@transitionMatrix, mcPrec@transitionMatrix)
  expect_identical(mc2@name, mcPrec@name)
})

test_that("the .yml extension is treated the same as .yaml", {
  skip_if_not_installed("yaml")
  tmp <- tempfile(fileext = ".yml")
  on.exit(unlink(tmp))

  toFile(mcRow, tmp)
  mc2 <- fromFile(tmp)
  expect_identical(mc2@transitionMatrix, mcRow@transitionMatrix)
})

## ---- toFile()/fromFile(): CSV round-trip and its documented limitation -

test_that("CSV round-trips the transition matrix but not the chain's name", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  toFile(mcRow, tmp)
  mc2 <- fromFile(tmp)

  expect_equal(unname(mc2@transitionMatrix), unname(mcRow@transitionMatrix))
  expect_equal(rownames(mc2@transitionMatrix), statesRow)
  expect_true(mc2@byrow)
  # No column in a CSV file can hold the chain's name, so it always comes
  # back as the package's own default rather than "RowChain".
  expect_equal(mc2@name, new("markovchain")@name)
})

test_that("fromFile() rejects a CSV whose header row and first column disagree", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  writeLines(c(',"a","b"', '"a",0.5,0.5', '"z",0.5,0.5'), tmp)
  expect_error(fromFile(tmp), "same state names")
})

## ---- round-trip agreement between toDictionary() and toFile() ----------

test_that("toFile()'s JSON/YAML output matches toDictionary() for the same chain", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("yaml")

  d <- toDictionary(mcRow)

  tmpJson <- tempfile(fileext = ".json")
  tmpYaml <- tempfile(fileext = ".yaml")
  on.exit(unlink(c(tmpJson, tmpYaml)))

  toFile(mcRow, tmpJson)
  toFile(mcRow, tmpYaml)

  expect_identical(fromFile(tmpJson)@transitionMatrix, fromDictionary(d)@transitionMatrix)
  expect_identical(fromFile(tmpYaml)@transitionMatrix, fromDictionary(d)@transitionMatrix)
})

## ---- missing optional dependencies (simulated) --------------------------

test_that("toFile()/fromFile() report an informative error when jsonlite/yaml is unavailable", {
  testthat::local_mocked_bindings(
    requireNamespace = function(...) FALSE,
    .package = "base"
  )
  expect_error(toFile(mcRow, tempfile(fileext = ".json")), "jsonlite")
  expect_error(toFile(mcRow, tempfile(fileext = ".yaml")), "yaml")
  expect_error(fromFile(tempfile(fileext = ".json")), "jsonlite|not found")
})
