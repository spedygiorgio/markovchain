context("selectOrder")

.extdata <- function(file, column) {
  path <- system.file("extdata", file, package = "markovchain")
  skip_if(!nzchar(path), paste(file, "not installed"))
  read.csv(path)[[column]]
}

# log-likelihood of a full order-k chain on observations t = start..n, by brute force
.bruteLogLik <- function(x, k, start) {
  x <- as.character(x)
  t <- start:length(x)
  ctx <- vapply(t, function(s) if (k == 0) "" else paste(x[(s - k):(s - 1)], collapse = "|"), "")
  ll <- 0
  for (c in unique(ctx)) {
    nxt <- x[t][ctx == c]
    freq <- table(nxt) / length(nxt)
    ll <- ll + sum(log(freq[nxt]))
  }
  ll
}

test_that("log-likelihoods, parameters and criteria follow their definitions", {
  set.seed(1)
  x <- sample(c("a", "b", "c"), 400, replace = TRUE)
  sel <- selectOrder(x, maxOrder = 3)
  tab <- sel$table
  expect_equal(tab$order, 0:3)
  expect_equal(sel$nobs, 397)
  for (k in 0:3) expect_equal(tab$logLik[k + 1], .bruteLogLik(x, k, 4))
  expect_equal(tab$npar, 3^(0:3) * 2)
  expect_equal(tab$AIC, -2 * tab$logLik + 2 * tab$npar)
  expect_equal(tab$BIC, -2 * tab$logLik + log(397) * tab$npar)
  expect_equal(tab$LR[-1], 2 * diff(tab$logLik))
  expect_equal(tab$df[-1], diff(tab$npar))
  expect_true(is.na(tab$LR[1]))
  expect_true(all(diff(tab$logLik) >= 0))
  # an independent sequence: BIC selects order 0
  expect_equal(sel$order, 0L)
})

test_that("the true order of a simulated chain is recovered", {
  set.seed(2)
  # order 2: the next state repeats the state two steps back with probability 0.85
  n <- 3000
  x <- character(n); x[1:2] <- c("a", "b")
  for (t in 3:n) x[t] <- if (runif(1) < 0.85) x[t - 2] else sample(c("a", "b", "c"), 1)
  expect_equal(selectOrder(x, maxOrder = 4)$order, 2L)
  expect_equal(selectOrder(x, maxOrder = 4, criterion = "AIC")$order, 2L)
  # first-order chain
  mc <- new("markovchain", states = c("a", "b", "c"),
            transitionMatrix = matrix(c(0.8, 0.1, 0.1, 0.2, 0.7, 0.1, 0.1, 0.2, 0.7), 3, byrow = TRUE))
  y <- rmarkovchain(3000, mc, t0 = "a")
  expect_equal(selectOrder(y, maxOrder = 3)$order, 1L)
})

test_that("selectOrder reproduces the Markov chain rows of Berchtold and Raftery (2002)", {
  wind <- .extdata("koeberg_wind.csv", "state")
  sel <- selectOrder(wind, maxOrder = 3, start = 15, parameters = "observed")
  expect_equal(sel$nobs, 730)
  expect_equal(sel$table$npar, c(3, 11, 27, 39))
  expect_lt(max(abs(sel$table$logLik - c(-954.8, -413.3, -374.9, -346.2))), 0.06)
  expect_lt(max(abs(sel$table$BIC - c(1929.4, 899.1, 927.9, 949.5))), 0.06)
  expect_equal(sel$order, 1L)
  seizures <- .extdata("epileptic_seizures.csv", "seizure")
  sel <- selectOrder(seizures, maxOrder = 3, start = 15, parameters = "observed")
  expect_equal(sel$table$npar, c(1, 2, 4, 8))
  expect_lt(max(abs(sel$table$logLik - c(-129.3, -122.6, -119.3, -117.2))), 0.06)
  expect_lt(max(abs(sel$table$BIC - c(263.9, 255.6, 259.6, 276.3))), 0.06)
  expect_equal(sel$order, 1L)
})

test_that("the likelihood-ratio test of order 1 against 0 equals assessIndependence", {
  wind <- as.character(.extdata("koeberg_wind.csv", "state"))
  sel <- selectOrder(wind, maxOrder = 1)
  g <- assessIndependence(wind, method = "G", verbose = FALSE)
  expect_equal(sel$table$LR[2], unname(g$statistic), tolerance = 1e-10)
  expect_equal(sel$table$df[2], unname(g$parameter))
})

test_that("a list of sequences pools the counts without crossing sequences", {
  set.seed(3)
  a <- sample(c("x", "y"), 50, TRUE); b <- sample(c("x", "y"), 70, TRUE)
  sel <- selectOrder(list(a, b), maxOrder = 2)
  expect_equal(sel$nobs, 48 + 68)
  # pooled counts: the log-likelihood of order 0 uses the pooled frequencies
  pooled <- c(a[3:50], b[3:70])
  expect_equal(sel$table$logLik[1], sum(log((table(pooled) / length(pooled))[pooled])))
  # order 1 contexts never cross from a to b
  t1 <- c(a[3:50], b[3:70]); c1 <- c(a[2:49], b[2:69])
  ll1 <- sum(vapply(unique(c1), function(cc) { v <- t1[c1 == cc]; sum(log((table(v) / length(v))[v])) }, 0))
  expect_equal(sel$table$logLik[2], ll1)
  # a single sequence in a list gives the same result as the sequence
  expect_equal(selectOrder(list(a), maxOrder = 2), selectOrder(a, maxOrder = 2))
})

test_that("numeric and factor sequences are accepted", {
  data(rain, package = "markovchain")
  ref <- selectOrder(rain$rain, maxOrder = 2)
  expect_equal(selectOrder(factor(rain$rain), maxOrder = 2), ref)
  expect_equal(selectOrder(as.integer(factor(rain$rain)), maxOrder = 2)$table, ref$table)
})

test_that("invalid input is rejected and sparse contexts are flagged", {
  x <- c("a", "b", "a", "b", "b", "a")
  expect_error(selectOrder(c("a", NA, "b")), "missing")
  expect_error(selectOrder(x, maxOrder = -1), "maxOrder")
  expect_error(selectOrder(x, maxOrder = 1.5), "maxOrder")
  expect_error(selectOrder(x, maxOrder = 2, start = 2), "start")
  expect_error(selectOrder(x, maxOrder = 2, start = 10), "no sequence")
  expect_error(selectOrder(rep("a", 10)), "two different states")
  expect_error(selectOrder(x, criterion = "foo"), "should be one of")
  expect_warning(selectOrder(sample(letters, 50, TRUE), maxOrder = 2), "unreliable")
})
