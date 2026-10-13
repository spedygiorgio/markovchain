context("Wood pewee song: published values and comparison with the march package")

.pewee <- function() {
  as.character(read.csv(system.file("extdata", "pewee_song.csv",
                                    package = "markovchain"))$phrase)
}

test_that("selectOrder reproduces the BIC of Kharin and Piatlitski (2011, Table 2)", {
  x <- .pewee()
  expect_length(x, 1327)
  # MC(s, s) for s = 1..4, all evaluated after the first 4 observations, zero
  # probabilities not counted as parameters
  sel <- selectOrder(x, maxOrder = 4, start = 5, parameters = "observed")
  expect_equal(round(sel$table$BIC[2:5], 1), c(1424.2, 801.9, 808.6, 768.3))
  expect_equal(sel$nobs, 1323)
})

test_that("fitMTD matches or improves on march::march.mtd.construct(maxOrder = 4)", {
  x <- .pewee()
  # log-likelihoods of march 3.3.2 (set.seed(1), default settings) on the same
  # 1323 observations: MTD(2) -568.023, MTD(3) -567.997, MTDg(2) -497.718,
  # MTDg(3) -490.673
  m2 <- fitMTD(x, 2, start = 5)
  m3 <- fitMTD(x, 3, start = 5)
  expect_equal(m2$logLikelihood, -567.996, tolerance = 1e-3, scale = 1)
  expect_equal(m3$logLikelihood, -567.996, tolerance = 1e-3, scale = 1)
  expect_equal(unname(m2$lambda), c(0, 1), tolerance = 1e-6)
  g2 <- fitMTD(x, 2, start = 5, model = "MTDg")
  g3 <- fitMTD(x, 3, start = 5, model = "MTDg")
  expect_gt(g2$logLikelihood, -497.718)
  expect_gt(g3$logLikelihood, -490.673)
  expect_equal(g2$logLikelihood, -486.155, tolerance = 1e-3, scale = 1)
  expect_equal(g3$logLikelihood, -483.985, tolerance = 1e-3, scale = 1)
  # the weight of lag 1 of MTD(2) is on the boundary
  expect_warning(fitMTD(x, 2, start = 5, se = TRUE), "boundary")
})
