test_that("Kullback et al. published example uses row-wise multinomial totals", {
  counts <- matrix(c(
    51, 11,  8,
    12, 31,  9,
     6, 11, 10
  ), nrow = 3, byrow = TRUE,
  dimnames = list(c("0", "1", "2"), c("0", "1", "2")))

  P <- matrix(c(
    5/8, 1/4, 1/8,
    1/4, 1/2, 1/4,
    1/4, 3/8, 3/8
  ), nrow = 3, byrow = TRUE,
  dimnames = dimnames(counts))

  mc <- new("markovchain", states = rownames(P), transitionMatrix = P)
  ans <- verifyEmpiricalToTheoretical(counts, mc, method = "G", verbose = FALSE)

  expect_s3_class(ans, "htest")
  expect_equal(unname(ans$statistic), 6.518384096, tolerance = 1e-8)
  expect_equal(unname(ans$parameter), 6)
  expect_equal(names(ans$statistic), "G-squared")
  expect_equal(names(ans$parameter), "df")
  expect_gt(ans$p.value, 0.05)
})
