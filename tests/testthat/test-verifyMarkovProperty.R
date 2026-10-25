test_that("verifyMarkovProperty returns the correct df for a binary chain", {
  sequence <- c(
    "a", "b", "a", "a", "a", "a", "b", "a", "b",
    "a", "b", "a", "a", "b", "b", "b", "a"
  )

  ans <- verifyMarkovProperty(
    sequence,
    method = "Pearson",
    verbose = FALSE
  )

  expect_equal(ans$dof, 2)
  expect_true(is.finite(ans$statistic))
  expect_true(is.finite(ans$p.value))
})


test_that("G and Pearson statistics are non-negative", {
  sequence <- c(
    "a", "b", "a", "a", "a", "a", "b", "a", "b",
    "a", "b", "a", "a", "b", "b", "b", "a"
  )

  g <- verifyMarkovProperty(sequence, method = "G", verbose = FALSE)
  p <- verifyMarkovProperty(sequence, method = "Pearson", verbose = FALSE)

  expect_gte(g$statistic, 0)
  expect_gte(p$statistic, 0)
  expect_equal(g$dof, p$dof)
})


test_that("simulation method is reproducible with a seed", {
  sequence <- c(
    "a", "b", "a", "a", "a", "a", "b", "a", "b",
    "a", "b", "a", "a", "b", "b", "b", "a"
  )

  ans1 <- verifyMarkovProperty(
    sequence,
    method = "simulation",
    B = 99,
    seed = 123,
    verbose = FALSE
  )
  ans2 <- verifyMarkovProperty(
    sequence,
    method = "simulation",
    B = 99,
    seed = 123,
    verbose = FALSE
  )

  expect_equal(ans1$p.value, ans2$p.value)
  expect_equal(ans1$simulations, ans2$simulations)
})


test_that("short sequences are rejected", {
  expect_error(
    verifyMarkovProperty(c("a", "b", "a"), verbose = FALSE),
    "at least four"
  )
})


test_that("missing values are rejected", {
  expect_error(
    verifyMarkovProperty(c("a", NA, "b", "a"), verbose = FALSE),
    "missing"
  )
})


test_that("simulation requires an outgoing transition for every state", {
  sequence <- c("a", "b", "b", "b", "c")

  expect_error(
    verifyMarkovProperty(
      sequence,
      method = "simulation",
      B = 9,
      verbose = FALSE
    ),
    "outgoing transition"
  )
})

test_that("the C++ simulation reproduces the former R loop exactly", {
  # The R code that method = "simulation" used until version 1.2: one call
  # to sample.int() per step and the G statistic computed in R.
  referenceSimulations <- function(sequence, B, seed) {
    states <- unique(sequence); r <- length(states); n <- length(sequence)
    idx <- match(sequence, states)
    triples <- function(x) array(tabulate(x[-c(1L, 2L)] + (x[2:(n - 1L)] - 1L) * r +
                                          (x[1:(n - 2L)] - 1L) * r * r, nbins = r^3L),
                                 dim = c(r, r, r))
    gStat <- function(N) {
      Nij <- apply(N, c(1L, 2L), sum); Njk <- apply(N, c(2L, 3L), sum); Nj <- rowSums(Njk)
      E <- array(0, dim = c(r, r, r))
      for (j in seq_len(r)) if (Nj[j] > 0L) E[, j, ] <- outer(Nij[, j], Njk[j, ]) / Nj[j]
      pos <- E > 0 & N > 0
      2 * sum(N[pos] * log(N[pos] / E[pos]))
    }
    N <- triples(idx)
    Njk <- apply(N, c(2L, 3L), sum)
    P <- Njk / rowSums(Njk)
    set.seed(seed)
    vapply(seq_len(B), function(b) {
      x <- integer(n); x[1L] <- idx[1L]
      for (t in seq_len(n - 1L)) x[t + 1L] <- sample.int(r, 1L, prob = P[x[t], ])
      gStat(triples(x))
    }, numeric(1))
  }
  # The sequences are the same everywhere (same random numbers, same
  # inversion). The statistics are bit for bit those of R on x86-64, where the
  # C++ sums are accumulated in the same 80-bit long double as R's sum(), and
  # agree to the last bits elsewhere (on macOS arm64 they differ in the last
  # bit, as seen in CI).
  sameBits <- R.version$arch %in% c("x86_64", "amd64")
  expectSame <- function(object, expected) {
    if (sameBits) expect_identical(object, expected) else expect_equal(object, expected, tolerance = 1e-12)
  }
  s <- c("a", "b", "c", "d")
  mc <- new("markovchain", states = s, transitionMatrix = matrix(
    c(0.5, 0.3, 0.2, 0, 0.1, 0.6, 0.3, 0, 0.25, 0.25, 0.25, 0.25, 0.2, 0.2, 0.3, 0.3),
    4, byrow = TRUE, dimnames = list(s, s)))
  for (seed in 1:3) {
    set.seed(seed)
    x <- rmarkovchain(150, mc, t0 = "a")
    ans <- verifyMarkovProperty(x, method = "simulation", B = 60, seed = seed, verbose = FALSE)
    expectSame(ans$simulations, referenceSimulations(x, 60, seed))
    expect_equal(ans$p.value, (1 + sum(ans$simulations >= ans$statistic)) / 61)
  }
  # the observed statistic, recomputed in C++ for the comparison, is the one of
  # method = "G"
  expectSame(ans$statistic, verifyMarkovProperty(x, method = "G", verbose = FALSE)$statistic)
  # ties between probabilities (a uniform row) and two states
  set.seed(4)
  y <- sample(c("u", "v"), 80, replace = TRUE)
  expectSame(verifyMarkovProperty(y, method = "simulation", B = 40, seed = 9,
                                  verbose = FALSE)$simulations,
             referenceSimulations(y, 40, 9))
})

test_that("a simulated sequence with the observed triple counts counts as a tie", {
  # In a short sequence of two states many simulated sequences have the same
  # triple counts as the observed one: their statistic must equal the observed
  # value exactly, on every platform, so that they enter the p-value (>=).
  x <- c("a", "b", "a", "a", "b", "b", "a", "b")
  ans <- verifyMarkovProperty(x, method = "simulation", B = 500, seed = 1, verbose = FALSE)
  expect_gt(sum(ans$simulations == ans$statistic), 0)
  expect_identical(unname(ans$statistic),
                   markovchain:::.markovPropertyStatisticRcpp(match(x, unique(x)), 2L, FALSE))
  expect_error(markovchain:::.markovPropertyStatisticRcpp(c(1L, 3L, 1L), 2L, FALSE), "codes")
})
