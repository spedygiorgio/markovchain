context("higherOrderPredict and higherOrderSimulate")

.wind <- function() {
  path <- system.file("extdata", "koeberg_wind.csv", package = "markovchain")
  skip_if(!nzchar(path), "koeberg_wind.csv not installed")
  read.csv(path)$state
}

test_that("predictions follow the definition of the MTD model", {
  x <- .wind()
  fit <- fitMTD(x, order = 2)
  Q <- fit$estimate@transitionMatrix
  p <- higherOrderPredict(fit, c(1, 2))
  # history oldest first: x_{t-2} = 1, x_{t-1} = 2
  expect_equal(p, fit$lambda[["lag1"]] * Q["2", ] + fit$lambda[["lag2"]] * Q["1", ])
  expect_equal(sum(p), 1)
  # only the last `order` states matter
  expect_equal(higherOrderPredict(fit, c(3, 4, 1, 2)), p)
})

test_that("predictions reproduce higherOrderLogLik", {
  x <- as.character(.wind())
  for (fit in list(fitMTD(x, order = 3), fitHigherOrder(x, order = 3, method = "mle"))) {
    times <- 4:length(x)
    histories <- t(vapply(times, function(t) x[(t - 3):(t - 1)], character(3)))
    probs <- higherOrderPredict(fit, histories)
    expect_equal(dim(probs), c(length(times), 4L))
    ll <- sum(log(probs[cbind(seq_along(times), match(x[times], colnames(probs)))]))
    expect_equal(ll, higherOrderLogLik(x, fit, start = 4)$logLik, tolerance = 1e-10)
  }
})

test_that("simulated sequences follow the predicted probabilities", {
  data(rain, package = "markovchain")
  fit <- fitHigherOrder(rain$rain, order = 2, method = "mle")
  set.seed(42)
  sim <- higherOrderSimulate(20000, fit, t0 = c("0", "0"), include.t0 = TRUE)
  expect_length(sim, 20002)
  expect_true(all(sim %in% rownames(fit$Q[[1]])))
  # empirical next-state frequencies after the history ("6+", "0")
  t <- which(sim[-c(length(sim) - 1, length(sim))] == "6+" &
             sim[-c(1, length(sim))] == "0") + 2
  freq <- prop.table(table(factor(sim[t], levels = names(fit$X))))
  expect_lt(max(abs(as.numeric(freq) - higherOrderPredict(fit, c("6+", "0")))), 0.03)
  expect_length(higherOrderSimulate(5, fit, t0 = c("0", "1-5")), 5)
  expect_length(higherOrderSimulate(0, fit, t0 = c("0", "1-5")), 0)
})

test_that("higherOrderSimulate is reproducible and deterministic chains are followed", {
  # a deterministic MTD(2): from a always b, from b always a, weight on lag 2 only
  Q <- matrix(c(0, 1, 1, 0), 2, dimnames = list(c("a", "b"), c("a", "b")))
  fit <- list(lambda = c(0, 1), Q = list(Q, Q))
  expect_equal(higherOrderSimulate(4, fit, t0 = c("a", "a")), c("b", "b", "a", "a"))
  set.seed(1); s1 <- higherOrderSimulate(10, fitMTD(.wind(), 2), t0 = c(1, 1))
  set.seed(1); s2 <- higherOrderSimulate(10, fitMTD(.wind(), 2), t0 = c(1, 1))
  expect_identical(s1, s2)
})

test_that("invalid input is rejected", {
  fit <- fitMTD(.wind(), order = 2)
  expect_error(higherOrderPredict(fit, 1), "at least 2 states")
  expect_error(higherOrderPredict(fit, c(1, 9)), "unknown")
  expect_error(higherOrderPredict(fit, c(1, NA)), "missing")
  expect_error(higherOrderPredict(list(a = 1), c(1, 2)), "fitHigherOrder")
  expect_error(higherOrderSimulate(-1, fit, t0 = c(1, 2)), "non-negative")
  expect_error(higherOrderSimulate(3, fit, t0 = 1), "at least 2 states")
})
