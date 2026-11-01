context("bootstrapProperties()")

data(rain, package = "markovchain")
rainChain <- new("markovchain", states = c("0", "1-5", "6+"),
                 transitionMatrix = createSequenceMatrix(rain$rain, toRowProbs = TRUE))

test_that("the estimates are the statistics of the maximum likelihood chain", {
  b <- bootstrapProperties(rain$rain, nboot = 20, seed = 1)
  expect_s3_class(b, "markovchainBootstrap")
  expect_equal(b$chain@transitionMatrix, rainChain@transitionMatrix)
  expect_equal(unname(b$estimate), c(slem(rainChain), spectralGap(rainChain),
                                     relaxationTime(rainChain), mixingTime(rainChain)))
  expect_identical(names(b$estimate), c("slem", "spectralGap", "relaxationTime", "mixingTime"))
  expect_equal(dim(b$replicates), c(20L, 4L))
  expect_equal(dim(b$ci), c(4L, 2L))
  expect_true(all(b$ci[, 1] <= b$ci[, 2]))
  expect_equal(b$counts, createSequenceMatrix(rain$rain))
  expect_output(print(b), "Percentile 95% intervals")
})

test_that("the multinomial bootstrap keeps the departures of each state", {
  n <- rowSums(createSequenceMatrix(rain$rain))
  b <- bootstrapProperties(rain$rain, nboot = 50, seed = 3, statistics = list(
    counts = function(mc) as.numeric(mc@transitionMatrix * n)))
  expect_equal(b$replicates, round(b$replicates), tolerance = 1e-9)
  rowTotals <- sapply(1:3, function(i) rowSums(b$replicates[, c(i, i + 3, i + 6)]))
  expect_equal(unname(rowTotals), matrix(n, 50, 3, byrow = TRUE))
})

test_that("the bootstrap standard error of a transition probability is binomial", {
  counts <- createSequenceMatrix(rain$rain)
  p <- counts[1, 2] / sum(counts[1, ])
  b <- bootstrapProperties(rain$rain, nboot = 2000, seed = 4, statistics = list(
    p12 = function(mc) mc@transitionMatrix[1, 2]))
  expect_equal(unname(b$se), sqrt(p * (1 - p) / sum(counts[1, ])), tolerance = 0.05)
  expect_lt(abs(unname(b$bias)), 0.003)
  # the parametric bootstrap, where the departures vary, gives a similar value
  bp <- bootstrapProperties(rain$rain, nboot = 500, seed = 4, method = "parametric",
                            statistics = list(p12 = function(mc) mc@transitionMatrix[1, 2]))
  expect_equal(unname(bp$se), unname(b$se), tolerance = 0.15)
})

test_that("seeds, interval types and input forms are consistent", {
  set.seed(99)
  before <- .Random.seed
  b1 <- bootstrapProperties(rain$rain, nboot = 30, seed = 5)
  expect_identical(.Random.seed, before)
  b2 <- bootstrapProperties(rain$rain, nboot = 30, seed = 5)
  expect_identical(b1$replicates, b2$replicates)
  basic <- bootstrapProperties(rain$rain, nboot = 30, seed = 5, type = "basic")
  expect_equal(unname(basic$ci), unname(cbind(2 * b1$estimate - b1$ci[, 2],
                                              2 * b1$estimate - b1$ci[, 1])))
  expect_equal(unname(bootstrapProperties(rain$rain, nboot = 30, seed = 5, level = 0.8)$ci),
               unname(t(apply(b1$replicates, 2, quantile, c(0.1, 0.9)))))
  # a count matrix gives the same multinomial bootstrap as the sequence
  bm <- bootstrapProperties(createSequenceMatrix(rain$rain), nboot = 30, seed = 5)
  expect_identical(bm$replicates, b1$replicates)
  # a list of sequences pools their transitions
  halves <- list(rain$rain[1:500], rain$rain[501:1096])
  bl <- bootstrapProperties(halves, nboot = 10, seed = 5)
  expect_equal(bl$counts, createSequenceMatrix(halves[[1]]) + createSequenceMatrix(halves[[2]]))
  expect_equal(nrow(bootstrapProperties(halves, nboot = 10, seed = 5, method = "parametric")$replicates),
               10L)
})

test_that("vector-valued and user statistics are flattened with names", {
  b <- bootstrapProperties(rain$rain, nboot = 20, seed = 6,
                           statistics = list("slem", dry = function(mc) steadyStates(mc)[1, "0"],
                                             pi = "steadyStates"))
  expect_identical(names(b$estimate), c("slem", "dry", "pi.0", "pi.1-5", "pi.6+"))
  expect_equal(sum(b$estimate[3:5]), 1)
  expect_equal(unname(rowSums(b$replicates[, 3:5])), rep(1, 20))
})

test_that("replicates where a statistic fails are counted and left out", {
  # state c is entered only once: many replicates make it unreachable, so the
  # chain is reducible and the spectral statistics cannot be computed
  x <- c(rep(c("a", "b"), 30), "c", "a", rep(c("b", "a"), 10))
  b <- bootstrapProperties(x, nboot = 100, seed = 7, statistics = list(
    "slem", "entropyRate", pAA = function(mc) mc@transitionMatrix["a", "a"]))
  expect_gt(b$failures[["slem"]], 0)
  # entropyRate() needs a unique stationary distribution too
  expect_equal(b$failures[["entropyRate"]], b$failures[["slem"]])
  expect_equal(b$failures[["pAA"]], 0)
  expect_true(all(is.finite(b$ci)))
  expect_output(print(b), "could not be computed")
})

test_that("invalid input is rejected", {
  expect_error(bootstrapProperties(rain$rain, nboot = 1), "nboot")
  expect_error(bootstrapProperties(rain$rain, level = 1), "level")
  expect_error(bootstrapProperties(rain$rain, statistics = "foo"), "unknown")
  expect_error(bootstrapProperties(rain$rain, statistics = list(function(mc) 1)), "names")
  expect_error(bootstrapProperties(rain$rain, statistics = list(a = 3)), "function")
  expect_error(bootstrapProperties(createSequenceMatrix(rain$rain), method = "parametric"),
               "sequences")
  expect_error(bootstrapProperties(c("a", "b", "c")), "never left")
  expect_error(bootstrapProperties(c("a", "b", NA, "a")), "missing")
  expect_error(bootstrapProperties(list("a")), "length at least 2")
  # a periodic estimated chain has no mixing time
  expect_error(bootstrapProperties(rep(c("a", "b"), 20)), "cannot be computed")
})

test_that("the C++ statistics equal slem(), spectralGap(), relaxationTime() and mixingTime()", {
  rstat <- function(mc) vapply(list(slem, spectralGap, relaxationTime, mixingTime), function(f)
    tryCatch(as.numeric(suppressWarnings(f(mc))), error = function(e) NA_real_), numeric(1))
  gen <- function(type, n) {
    M <- switch(type,
      dense = matrix(rexp(n * n), n),
      sparse = { M <- matrix(rexp(n * n) * (runif(n * n) < 0.35), n)
                 diag(M) <- diag(M) + (rowSums(M) == 0); M },
      cycle = { M <- matrix(0, n, n); M[cbind(1:n, c(2:n, 1))] <- 1; M },          # periodic
      reducible = { M <- matrix(rexp(n * n), n); M[1, ] <- 0; M[1, 1] <- 1; M },
      symmetric = { M <- matrix(rexp(n * n), n); M + t(M) },
      rotation = { M <- matrix(0, n, n); M[cbind(1:n, c(2:n, 1))] <- 0.8           # complex
                   diag(M) <- 0.2; M },                                             # eigenvalues
      lazy = { M <- matrix(rexp(n * n), n); 0.97 * diag(n) + 0.03 * M / rowSums(M) })
    M <- M / rowSums(M)
    s <- paste0("s", seq_len(n))
    dimnames(M) <- list(s, s)
    new("markovchain", states = s, transitionMatrix = M)
  }
  set.seed(2026)
  for (type in c("dense", "sparse", "cycle", "reducible", "symmetric", "rotation", "lazy")) {
    for (n in c(1, 2, 3, 6)) {
      if (n == 1 && type %in% c("cycle", "rotation")) next
      for (r in 1:3) {
        mc <- gen(type, n)
        expect_equal(markovchain:::.spectralStatisticsRcpp(mc@transitionMatrix, 0:3, 0.25, 10000L),
                     rstat(mc), tolerance = 1e-12, info = paste(type, n))
      }
    }
  }
  # the stopping rules of mixingTime(): epsilon and maxIter
  lazy <- gen("lazy", 4)
  expect_equal(markovchain:::.spectralStatisticsRcpp(lazy@transitionMatrix, 3L, 0.01, 10000L),
               mixingTime(lazy, epsilon = 0.01))
  expect_true(is.na(markovchain:::.spectralStatisticsRcpp(lazy@transitionMatrix, 3L, 0.01, 2L)))
})

test_that("the C++ and R paths give the same replicates", {
  wrapped <- list(slem = function(m) slem(m), spectralGap = function(m) spectralGap(m),
                  relaxationTime = function(m) relaxationTime(m),
                  mixingTime = function(m) mixingTime(m))
  for (method in c("multinomial", "parametric")) {
    fast <- bootstrapProperties(rain$rain, nboot = 100, method = method, seed = 11)
    slow <- bootstrapProperties(rain$rain, nboot = 100, method = method, seed = 11,
                                statistics = wrapped)
    expect_equal(fast$replicates, slow$replicates, tolerance = 1e-12, info = method)
  }
  # a built-in under another name still uses its definition
  renamed <- bootstrapProperties(rain$rain, nboot = 20, seed = 3,
                                 statistics = list(gap = "spectralGap"))
  expect_equal(unname(renamed$replicates[, 1]),
               unname(bootstrapProperties(rain$rain, nboot = 20, seed = 3,
                                          statistics = "spectralGap")$replicates[, 1]))
})

test_that("the multinomial counts are those of stats::rmultinom(), state by state", {
  counts <- createSequenceMatrix(rain$rain)
  P <- counts / rowSums(counts)
  set.seed(21)
  expected <- vapply(1:25, function(b) {
    N <- t(vapply(1:3, function(i) as.numeric(stats::rmultinom(1L, sum(counts[i, ]), P[i, ])),
                  numeric(3)))
    mc <- new("markovchain", states = rownames(P), transitionMatrix = N / rowSums(N))
    slem(mc)
  }, numeric(1))
  b <- bootstrapProperties(rain$rain, nboot = 25, seed = 21, statistics = "slem")
  expect_equal(unname(b$replicates[, 1]), expected, tolerance = 1e-12)
})
