context("Missing optional packages are signalled with a warning")

test_that("the msm coercion, and the diagram and DiagrammeR helpers, warn instead of printing", {
  testthat::local_mocked_bindings(requireNamespace = function(package, ..., quietly = FALSE) FALSE,
                                  .package = "base")
  mc <- new("markovchain", states = c("a", "b"), transitionMatrix = matrix(c(.5, .5, .2, .8), 2, byrow = TRUE))
  expect_warning(out <- markovchain:::.msm2Mc(NULL), "msm")
  expect_null(out)
  expect_silent(suppressWarnings(markovchain:::.msm2Mc(NULL)))
  expect_warning(out <- markovchain:::.plotdiagram(mc), "diagram")
  expect_null(out)
  expect_warning(out <- markovchain:::.plotDiagrammeR(mc), "DiagrammeR")
  expect_null(out)
  # nothing is written to the output stream any more
  expect_output(suppressWarnings(markovchain:::.msm2Mc(NULL)), NA)
})
