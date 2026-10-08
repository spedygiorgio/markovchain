context("XML in toFile()/fromFile() and summary(details = TRUE)")

s <- c("sun", "rain")
mc <- new("markovchain", states = s, name = "Weather <&> \"test\"",
          transitionMatrix = matrix(c(0.7, 0.3, 0.45, 0.55), 2, byrow = TRUE,
                                    dimnames = list(s, s)))

test_that("XML round trip is exact and keeps the name", {
  skip_if_not_installed("xml2")
  f <- tempfile(fileext = ".xml")
  toFile(mc, f)
  back <- fromFile(f)
  expect_identical(back@transitionMatrix, mc@transitionMatrix)
  expect_identical(back@states, s)
  expect_identical(back@name, mc@name)
  # irrational probabilities survive too
  r <- randomMarkovChain(5, seed = 8)
  toFile(r, f)
  expect_identical(fromFile(f)@transitionMatrix, r@transitionMatrix)
  # byrow = FALSE chains are written by rows
  cc <- new("markovchain", states = s, transitionMatrix = t(mc@transitionMatrix),
            byrow = FALSE)
  toFile(cc, f, format = "xml")
  expect_true(fromFile(f)@byrow)
  expect_equal(fromFile(f)@transitionMatrix, mc@transitionMatrix)
  unlink(f)
})

test_that("a file written by PyDTMC is read", {
  skip_if_not_installed("xml2")
  f <- tempfile(fileext = ".xml")
  writeLines(c("<?xml version='1.0' encoding='utf-8' standalone='yes' ?>",
               "<MarkovChain>",
               "\t<Item state_from=\"sun\" state_to=\"sun\" probability=\"0.7\"/>",
               "\t<Item state_from=\"sun\" state_to=\"rain\" probability=\"0.3\"/>",
               "\t<Item state_from=\"rain\" state_to=\"sun\" probability=\"0.45\"/>",
               "\t<Item state_from=\"rain\" state_to=\"rain\" probability=\"0.55\"/>",
               "</MarkovChain>"), f)
  py <- fromFile(f)
  expect_identical(py@states, s)
  expect_equal(py@transitionMatrix, mc@transitionMatrix)
  expect_identical(py@name, "Unnamed Markov chain")
  unlink(f)
})

test_that("malformed XML files are rejected", {
  skip_if_not_installed("xml2")
  f <- tempfile(fileext = ".xml")
  writeLines("<Chain><Item state_from=\"a\" state_to=\"a\" probability=\"1\"/></Chain>", f)
  expect_error(fromFile(f), "MarkovChain")
  writeLines("<MarkovChain><Item state_from=\"a\" state_to=\"a\"/></MarkovChain>", f)
  expect_error(fromFile(f), "attributes")
  writeLines(c("<MarkovChain>",
               "<Item state_from=\"a\" state_to=\"a\" probability=\"0.5\"/>",
               "<Item state_from=\"a\" state_to=\"b\" probability=\"0.5\"/>",
               "<Item state_from=\"b\" state_to=\"b\" probability=\"1\"/>",
               "</MarkovChain>"), f)
  expect_error(fromFile(f), "every pair")
  writeLines("<MarkovChain><Item state_from=\"a\" state_to=\"a\" probability=\"x\"/></MarkovChain>", f)
  expect_error(fromFile(f), "not numbers")
  unlink(f)
})

test_that("summary() without details is unchanged; details adds properties", {
  plain <- capture.output(res <- summary(mc))
  expect_null(res$details)
  detailed <- capture.output(resd <- summary(mc, details = TRUE))
  expect_identical(detailed[seq_along(plain)], plain)
  expect_true(any(grepl("Further properties", detailed)))
  d <- resd$details
  expect_equal(d$size, 2L)
  expect_true(d$irreducible)
  expect_equal(d$period, 1)
  expect_true(d$reversible)        # every 2-state irreducible chain is
  expect_false(d$absorbingChain)
  expect_equal(d$entropyRate, entropyRate(mc))
  expect_equal(d$kemenyConstant, kemenyConstant(mc))
  expect_error(summary(mc, details = NA), "details")
})

test_that("summary details on reducible and absorbing chains", {
  gr <- gamblersRuin(4, 0.5)
  d <- markovchain:::.summaryDetails(gr)
  expect_true(d$absorbingChain)
  expect_false(d$irreducible)
  expect_true(is.na(d$period))
  expect_true(is.na(d$kemenyConstant))
  expect_true(is.na(d$entropyRate))  # two recurrent classes
  id <- markovchain:::.summaryDetails(identityChain(3))
  expect_true(id$symmetric)
  expect_equal(id$classes, 3L)
  expect_equal(id$rank, 3L)
})
