test_that("autoplot.markovchain returns a ggplot object", {
  skip_if_not_installed("ggplot2")

  P <- matrix(
    c(0.7, 0.2, 0.1,
      0.3, 0.4, 0.3,
      0.2, 0.45, 0.35),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("sunny", "cloudy", "rain"),
                    c("sunny", "cloudy", "rain"))
  )
  mc <- new("markovchain", transitionMatrix = P, name = "Weather")

  plt <- ggplot2::autoplot(mc)

  expect_s3_class(plt, "ggplot")
  expect_equal(plt$labels$title, "Weather")
})

test_that("autoplot.markovchain handles structural zeros and suppresses labels", {
  skip_if_not_installed("ggplot2")

  P <- matrix(
    c(1, 0, 0,
      0, 0.5, 0.5,
      0, 1, 0),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(letters[1:3], letters[1:3])
  )
  mc <- new("markovchain", transitionMatrix = P)

  plt <- ggplot2::autoplot(mc, show_probabilities = FALSE)

  expect_s3_class(plt, "ggplot")
  expect_false(any(vapply(plt$layers, function(layer) {
    inherits(layer$geom, "GeomLabel")
  }, logical(1))))
})
