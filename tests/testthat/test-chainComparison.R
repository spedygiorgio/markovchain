context("compareChains() and divergenceRate()")

s <- c("a", "b", "c")
P <- new("markovchain", states = s, transitionMatrix = matrix(
  c(0.5, 0.3, 0.2, 0.1, 0.6, 0.3, 0.25, 0.25, 0.5), 3, byrow = TRUE,
  dimnames = list(s, s)))
Q <- new("markovchain", states = s, transitionMatrix = matrix(
  c(0.4, 0.4, 0.2, 0.2, 0.5, 0.3, 0.3, 0.1, 0.6), 3, byrow = TRUE,
  dimnames = list(s, s)))

test_that("the measures agree with PyDTMC 10.0.0", {
  # calculate_transition_distance(), calculate_transition_divergence() and
  # calculate_stationary_distance() on the same two chains
  cc <- compareChains(P, Q)
  expect_equal(cc$distances[["frobenius"]], 0.27386127875258304, tolerance = 1e-14)
  expect_equal(cc$distances[["hellinger"]], 0.14248607206963021, tolerance = 1e-14)
  expect_equal(cc$distances[["maximum"]], 0.15, tolerance = 1e-14)
  expect_equal(cc$distances[["totalVariation"]], 0.15, tolerance = 1e-14)
  expect_equal(cc$distances[["stationary"]], 0.088250930356193608, tolerance = 1e-12)
  expect_equal(cc$divergences[["maxRowXY"]], 0.092331515373072839, tolerance = 1e-14)
  expect_equal(cc$divergences[["jensenShannon"]], 0.020017380814281917, tolerance = 1e-14)
  expect_null(cc$test)
  expect_s3_class(cc, "markovchainComparison")
  expect_output(print(cc), "Distances")
})

test_that("divergenceRate() matches the closed form for two states", {
  two <- function(a, b) new("markovchain", states = c("0", "1"), transitionMatrix =
    matrix(c(1 - a, a, b, 1 - b), 2, byrow = TRUE, dimnames = list(c("0", "1"), c("0", "1"))))
  kl <- function(p, q) p * log(p / q) + (1 - p) * log((1 - p) / (1 - q))
  a <- 0.3; b <- 0.1; c <- 0.5; d <- 0.2
  pi0 <- b / (a + b)
  expected <- pi0 * kl(a, c) + (1 - pi0) * kl(b, d)
  expect_equal(divergenceRate(two(a, b), two(c, d)), expected, tolerance = 1e-14)
  expect_equal(divergenceRate(two(a, b), two(c, d), base = 2), expected / log(2),
               tolerance = 1e-14)
  expect_equal(divergenceRate(P, P), 0)
  expect_gt(divergenceRate(P, Q), 0)
  expect_false(isTRUE(all.equal(divergenceRate(P, Q), divergenceRate(Q, P))))
})

test_that("divergenceRate() is the limit of the log-likelihood ratio per step", {
  set.seed(2026)
  x <- rmarkovchain(1e5, P, t0 = "a")
  i <- match(x[-length(x)], s)
  j <- match(x[-1], s)
  llr <- mean(log(P@transitionMatrix[cbind(i, j)] / Q@transitionMatrix[cbind(i, j)]))
  expect_equal(llr, divergenceRate(P, Q), tolerance = 0.03)
})

test_that("2 N times the rate with empirical weights is Kullback's statistic", {
  # Kullback, Kupperman and Ku (1962), Table 9.2: homogeneity of two chains,
  # 2I = 160.961 (also computed by verifyHomogeneity()).
  data(kullback, package = "markovchain")
  pooled <- kullback[[1]] + kullback[[2]]
  pbar <- new("markovchain", states = rownames(pooled),
              transitionMatrix = pooled / rowSums(pooled))
  stat <- sum(vapply(kullback, function(m) {
    fit <- new("markovchain", states = rownames(m), transitionMatrix = m / rowSums(m))
    2 * sum(m) * divergenceRate(fit, pbar, weights = rowSums(m))
  }, numeric(1)))
  expect_equal(round(stat, 3), 160.961)
  expect_equal(stat, unname(verifyHomogeneity(kullback, verbose = FALSE)$statistic),
               tolerance = 1e-12)
  # and G^2 of verifyEmpiricalToTheoretical() for a fully specified chain
  counts <- matrix(c(51, 11, 8, 12, 31, 9, 6, 11, 10), 3, byrow = TRUE,
                   dimnames = list(s, s))
  theory <- new("markovchain", states = s, transitionMatrix = matrix(
    c(5/8, 1/4, 1/8, 1/4, 1/2, 1/4, 1/4, 3/8, 3/8), 3, byrow = TRUE, dimnames = list(s, s)))
  fit <- new("markovchain", states = s, transitionMatrix = counts / rowSums(counts))
  g2 <- 2 * sum(counts) * divergenceRate(fit, theory, weights = rowSums(counts))
  expect_equal(g2, unname(verifyEmpiricalToTheoretical(counts, theory, verbose = FALSE)$statistic),
               tolerance = 1e-12)
})

test_that("states are matched by name and byrow = FALSE is handled", {
  perm <- c("c", "a", "b")
  Qperm <- new("markovchain", states = perm,
               transitionMatrix = Q@transitionMatrix[perm, perm])
  expect_equal(compareChains(P, Qperm)$distances, compareChains(P, Q)$distances)
  Qcol <- new("markovchain", states = s, byrow = FALSE, transitionMatrix = t(Q@transitionMatrix))
  expect_equal(divergenceRate(P, Qcol), divergenceRate(P, Q))
  expect_equal(compareChains(P, Qcol)$divergences, compareChains(P, Q)$divergences)
})

test_that("infinite divergences, zero weights and several stationary distributions", {
  # Q forbids a -> c, which P allows
  Z <- new("markovchain", states = s, transitionMatrix = matrix(
    c(0.5, 0.5, 0, 0.1, 0.6, 0.3, 0.25, 0.25, 0.5), 3, byrow = TRUE, dimnames = list(s, s)))
  expect_identical(divergenceRate(P, Z), Inf)
  expect_true(is.finite(divergenceRate(Z, P)))
  expect_true(is.finite(compareChains(P, Z)$divergences[["jensenShannon"]]))
  # a zero weight removes the infinite row
  expect_true(is.finite(divergenceRate(P, Z, weights = c(a = 0, b = 1, c = 1))))
  # reducible chain with two closed classes: no unique stationary distribution
  R <- new("markovchain", states = s, transitionMatrix = matrix(
    c(1, 0, 0, 0, 1, 0, 0.5, 0.25, 0.25), 3, byrow = TRUE, dimnames = list(s, s)))
  expect_error(divergenceRate(R, P), "stationary distributions")
  expect_true(is.finite(divergenceRate(R, P, weights = c(1, 1, 1))))
  cc <- compareChains(R, P)
  expect_true(is.na(cc$divergences[["rateXY"]]))
  expect_true(is.na(cc$distances[["stationary"]]))
  # P allows a -> b and a -> c, which R forbids
  expect_identical(cc$divergences[["rateYX"]], Inf)
})

test_that("two sequences are estimated on the same states and tested", {
  set.seed(1)
  x <- rmarkovchain(400, P, t0 = "a")
  y <- rmarkovchain(400, P, t0 = "b")
  z <- rmarkovchain(400, Q, t0 = "a")
  same <- compareChains(x, y)
  diff <- compareChains(x, z)
  expect_s3_class(same$test, "htest")
  expect_gt(same$test$p.value, 0.05)
  expect_lt(diff$test$p.value, 0.01)
  expect_equal(unname(same$test$statistic),
               unname(verifyHomogeneity(list(x, y), verbose = FALSE)$statistic))
  expect_equal(same$x@transitionMatrix,
               createSequenceMatrix(x, toRowProbs = TRUE, possibleStates = s)[s, s],
               check.attributes = FALSE)
  expect_null(compareChains(x, z, test = "none")$test)
  expect_output(print(diff), "same chain")
  # a state never left in one sequence
  expect_error(compareChains(c(x, "d"), y), "never left")
})

test_that("invalid input is rejected", {
  other <- new("markovchain", states = c("a", "b"), transitionMatrix = diag(2))
  expect_error(compareChains(P, other), "same states")
  expect_error(compareChains(P, c("a", "b")), "both")
  expect_error(divergenceRate(P, Q, base = 1), "base")
  expect_error(divergenceRate(P, Q, weights = c(1, -1, 1)), "weights")
  expect_error(divergenceRate(P, Q, weights = c(x = 1, y = 1, z = 1)), "names")
})
