context("plotSequence()")

s <- c("sunny", "cloudy", "rain")
weather <- new("markovchain", states = s, name = "weather", transitionMatrix = matrix(
  c(0.7, 0.2, 0.1, 0.3, 0.4, 0.3, 0.2, 0.45, 0.35), 3, byrow = TRUE,
  dimnames = list(s, s)))

# draw on a null device
withNullDevice <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expr
}

test_that("simulated sequences have the requested shape and are reproducible", {
  out <- withNullDevice(plotSequence(weather, nSequences = 7, steps = 12, initial = "rain",
                                     seed = 1))
  expect_equal(dim(out$sequences), c(7L, 12L))
  expect_true(all(out$sequences[, 1] == "rain"))
  expect_true(all(out$sequences %in% s))
  set.seed(5)
  before <- .Random.seed
  again <- withNullDevice(plotSequence(weather, nSequences = 7, steps = 12, initial = "rain",
                                       seed = 1))
  expect_identical(.Random.seed, before)
  expect_identical(again$sequences, out$sequences)
  # every observed transition is possible in the chain
  from <- match(out$sequences[, -12], s)
  to <- match(out$sequences[, -1], s)
  expect_true(all(weather@transitionMatrix[cbind(from, to)] > 0))
})

test_that("proportions are the observed frequencies and expected follows redistribute()", {
  out <- withNullDevice(plotSequence(weather, type = "distribution", nSequences = 50,
                                     steps = 6, initial = c(0.5, 0.5, 0), seed = 2))
  manual <- sapply(1:6, function(t) table(factor(out$sequences[, t], levels = s)) / 50)
  expect_equal(unname(out$proportions), unname(manual))
  expect_identical(rownames(out$proportions), s)
  expect_equal(unname(t(out$expected)),
               unname(redistribute(weather, steps = 5, initial = c(0.5, 0.5, 0))))
  # law of large numbers on the last time point
  big <- withNullDevice(plotSequence(weather, type = "distribution", nSequences = 4000,
                                     steps = 4, seed = 3))
  expect_lt(max(abs(big$proportions[, 4] - big$expected[, 4])), 0.03)
})

test_that("observed data of several shapes are accepted", {
  data(holson, package = "markovchain")
  out <- withNullDevice(plotSequence(holson[1:30, 2:12], sort = "first"))
  expect_equal(dim(out$sequences), c(30L, 11L))
  expect_false(is.unsorted(as.numeric(out$sequences[, 1])))
  expect_null(out$expected)
  last <- withNullDevice(plotSequence(holson[1:30, 2:12], sort = "last"))
  expect_false(is.unsorted(as.numeric(last$sequences[, 11])))
  # a list of sequences of different lengths is padded with NA
  ragged <- withNullDevice(plotSequence(list(c("a", "b", "a"), c("b", "b"), "a")))
  expect_equal(dim(ragged$sequences), c(3L, 3L))
  expect_true(is.na(ragged$sequences[3, 2]))
  expect_equal(unname(ragged$proportions[, 2]), c(0, 1))   # a, b at time 2
  expect_equal(unname(ragged$proportions[, 3]), c(1, 0))
  # factor levels give the order of the states, unobserved levels included
  f <- factor(c("lo", "hi", "hi"), levels = c("lo", "mid", "hi"))
  expect_identical(rownames(withNullDevice(plotSequence(f))$proportions), c("lo", "mid", "hi"))
  # a fitted chain gives the expected distribution from the observed start
  fitted <- markovchainFit(holson[, 2:12])$estimate
  d <- withNullDevice(plotSequence(holson[, 2:12], type = "distribution", chain = fitted))
  expect_equal(d$expected[, 1], d$proportions[, 1])
  expect_equal(unname(colSums(d$expected)), rep(1, 11))
})

test_that("graphical parameters are restored", {
  withNullDevice({
    before <- par("mar")
    plotSequence(weather, nSequences = 3, steps = 3, seed = 1)
    expect_equal(par("mar"), before)
  })
})

test_that("invalid input is rejected", {
  withNullDevice({
    expect_error(plotSequence(weather, nSequences = 0), "nSequences")
    expect_error(plotSequence(weather, initial = "snow"), "state")
    expect_error(plotSequence(list()), "markovchain")
    expect_error(plotSequence(c("a", "b"), chain = weather), "not in the chain")
    expect_error(plotSequence(c("a", "b"), chain = "x"), "markovchain")
    expect_error(plotSequence(weather, col = "red"), "colour")
  })
})
