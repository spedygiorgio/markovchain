test_that("sparse transition matrix is accepted", {
  expect_is(as(sparsematrix, "markovchain"), "markovchain")
})
