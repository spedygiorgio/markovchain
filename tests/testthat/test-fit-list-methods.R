library(testthat)
library(markovchain)

# Which fitting methods accept a list of sequences (issue #165).
#
# "mle" and "map" always did. "laplace" now does too: transition counts are
# pooled over the sequences by createSequenceMatrix() and then smoothed, which
# is the same operation as for a single sequence. "bootstrap" is still refused,
# but with a message that says why rather than the old "method not available
# for a list".

c1 <- c("a", "b", "c", "c", "e")
c2 <- c("a", "b", "d", "e")
c3 <- c("a", "c", "b", "c", "d")
c4 <- c("a", "b", "b", "d", "b", "c", "d", "e")
c5 <- c("a", "c", "c", "d", "d")
c6 <- c("a", "c", "d", "d", "b", "b", "e")
mylist <- list(c1, c2, c3, c4, c5, c6)

test_that("laplace fitting accepts a list of sequences", {
  fit <- markovchainFit(data = mylist, method = "laplace", laplacian = 0.5)
  expect_s4_class(fit$estimate, "markovchain")
  expect_true(validObject(fit$estimate))
})

test_that("laplace on a list smooths the pooled transition counts", {
  laplacian <- 0.5
  fit <- markovchainFit(data = mylist, method = "laplace", laplacian = laplacian)

  counts <- createSequenceMatrix(mylist)
  expected <- (counts + laplacian) / rowSums(counts + laplacian)

  expect_equal(fit$estimate@transitionMatrix, expected,
               tolerance = 1e-12, check.attributes = FALSE)
})

test_that("laplace on a one-element list matches the same single sequence", {
  fromList <- markovchainFit(data = list(c1), method = "laplace",
                             laplacian = 0.3)$estimate
  fromSeq <- markovchainFit(data = c1, method = "laplace",
                            laplacian = 0.3)$estimate
  expect_equal(fromList@transitionMatrix, fromSeq@transitionMatrix,
               tolerance = 1e-12)
})

test_that("laplace on a list honours byrow", {
  byRow <- markovchainFit(data = mylist, method = "laplace", laplacian = 0.5,
                          byrow = TRUE)$estimate
  byCol <- markovchainFit(data = mylist, method = "laplace", laplacian = 0.5,
                          byrow = FALSE)$estimate

  expect_true(byRow@byrow)
  expect_false(byCol@byrow)
  expect_true(validObject(byCol))
  expect_equal(byCol@transitionMatrix, t(byRow@transitionMatrix),
               tolerance = 1e-12)
})

test_that("laplace with zero smoothing on a list reduces to the mle fit", {
  # With laplacian = 0 there is nothing to smooth, so the pooled counts are
  # just row-normalized -- which is exactly the maximum likelihood estimate.
  laplaceFit <- markovchainFit(data = mylist, method = "laplace",
                               laplacian = 0, sanitize = TRUE)$estimate
  mleFit <- markovchainFit(data = mylist, method = "mle",
                           sanitize = TRUE)$estimate
  expect_equal(laplaceFit@transitionMatrix, mleFit@transitionMatrix,
               tolerance = 1e-12)
})

test_that("bootstrap on a list is refused with an explanatory message", {
  expect_error(markovchainFit(data = mylist, method = "bootstrap", nboot = 5),
               "not supported for a list")
  # The message should point the user somewhere useful.
  expect_error(markovchainFit(data = mylist, method = "bootstrap", nboot = 5),
               "mle")
})

test_that("mle and map still accept lists", {
  for (method in c("mle", "map")) {
    fit <- markovchainFit(data = mylist, method = method, sanitize = TRUE)
    expect_s4_class(fit$estimate, "markovchain")
    expect_true(validObject(fit$estimate))
  }
})
