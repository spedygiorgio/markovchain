library(testthat)
library(markovchain)

randomChain <- function(n, seed) {
  set.seed(seed)
  M <- matrix(runif(n * n), n, n)
  M <- M / rowSums(M)
  nm <- letters[seq_len(n)]
  dimnames(M) <- list(nm, nm)
  new("markovchain", states = nm, transitionMatrix = M)
}

# A chain exactly lumpable into {a,b} / {c,d}: within each block, every state
# has the same total probability of moving to the other block (and to
# itself), so no information is lost by only tracking which block you are in.
statesLump <- c("a", "b", "c", "d")
Plump <- matrix(c(0.1, 0.1, 0.4, 0.4,
                   0.1, 0.1, 0.4, 0.4,
                   0.3, 0.3, 0.2, 0.2,
                   0.3, 0.3, 0.2, 0.2), byrow = TRUE, nrow = 4,
                 dimnames = list(statesLump, statesLump))
mcLump <- new("markovchain", states = statesLump, transitionMatrix = Plump)

## ---- exact lumpability is recovered with zero divergence -------------------

test_that("an exactly lumpable chain is aggregated with (near) zero KL divergence", {
  expect_true(is.lumpable(mcLump, list(g1 = c("a", "b"), g2 = c("c", "d"))))

  for (method in c("spectral-top-down", "spectral-bottom-up")) {
    result <- aggregateStates(mcLump, k = 2, method = method)
    expect_equal(result$klDivergence, 0, tolerance = 1e-8, info = method)

    grouped <- unname(sort(vapply(result$partition, function(s) paste(sort(s), collapse = ""), "")))
    expect_equal(grouped, sort(c("ab", "cd")), info = method)
    expect_true(validObject(result$aggregatedChain), info = method)
    expect_equal(unname(rowSums(result$aggregatedChain@transitionMatrix)),
                 c(1, 1), tolerance = 1e-10, info = method)
  }
})

test_that("the aggregated chain matches the exact lump() over the same partition", {
  result <- aggregateStates(mcLump, k = 2, method = "spectral-top-down")
  exact <- lump(mcLump, result$partition, force = FALSE)
  # Both use "Macro_1"/"Macro_2" states in the same order (partition list
  # order), so the transition matrices should agree up to relabeling.
  expect_equal(unname(result$aggregatedChain@transitionMatrix),
               unname(exact@transitionMatrix), tolerance = 1e-8)
})

test_that("strong lumpability alone does not imply zero KL divergence", {
  # The classical Land of Oz chain: {rainy, snowy} -> Bad_Weather is exactly
  # (Kemeny-Snell) lumpable -- both send a total of 0.75 into the block --
  # but rainy and snowy split that 0.75 differently between themselves
  # (0.5/0.25 vs 0.25/0.5), so aggregateStates()'s stricter KL divergence is
  # NOT zero here, even though it still recovers the same partition and the
  # same aggregated chain as lump(). This is a documented property of the
  # measure, not a bug: a zero divergence implies strong lumpability, but
  # not conversely.
  statesOz <- c("rainy", "nice", "snowy")
  Poz <- matrix(c(0.5, 0.25, 0.25,
                  0.5, 0.00, 0.5,
                  0.25, 0.25, 0.5), byrow = TRUE, nrow = 3,
                dimnames = list(statesOz, statesOz))
  mcOz <- new("markovchain", states = statesOz, transitionMatrix = Poz)
  partitionOz <- list(Bad_Weather = c("rainy", "snowy"), Nice_Weather = "nice")

  expect_true(is.lumpable(mcOz, partitionOz))

  result <- aggregateStates(mcOz, k = 2)
  grouped <- unname(sort(vapply(result$partition, function(s) paste(sort(s), collapse = ""), "")))
  expect_equal(grouped, sort(c("nice", "rainysnowy")))
  expect_gt(result$klDivergence, 1e-4)

  exact <- lump(mcOz, partitionOz, force = FALSE)
  expect_equal(unname(result$aggregatedChain@transitionMatrix),
               unname(exact@transitionMatrix), tolerance = 1e-8)
})

## ---- structural invariants on a generic chain -------------------------------

test_that("the returned partition covers every state exactly once, for both methods", {
  mc <- randomChain(10, seed = 1)
  for (method in c("spectral-top-down", "spectral-bottom-up")) {
    for (k in c(2, 5, 9)) {
      result <- aggregateStates(mc, k = k, method = method)
      allStates <- unname(unlist(result$partition))
      expect_equal(sort(allStates), sort(states(mc)), info = paste(method, k))
      expect_length(allStates, 10)
      expect_equal(length(result$partition), k, info = paste(method, k))
      expect_true(validObject(result$aggregatedChain), info = paste(method, k))
      expect_equal(unname(rowSums(result$aggregatedChain@transitionMatrix)),
                   rep(1, k), tolerance = 1e-8, info = paste(method, k))
    }
  }
})

test_that("the KL divergence is non-negative", {
  mc <- randomChain(9, seed = 5)
  for (method in c("spectral-top-down", "spectral-bottom-up")) {
    for (k in c(2, 4, 8)) {
      result <- aggregateStates(mc, k = k, method = method)
      expect_gte(result$klDivergence, -1e-8, label = paste(method, k))
    }
  }
})

test_that("column-stochastic storage gives the same divergence as row-stochastic", {
  mc <- randomChain(7, seed = 3)
  mcCol <- new("markovchain", states = states(mc), byrow = FALSE,
               transitionMatrix = t(mc@transitionMatrix))

  byRow <- aggregateStates(mc, k = 3, method = "spectral-bottom-up")
  byCol <- aggregateStates(mcCol, k = 3, method = "spectral-bottom-up")
  expect_equal(byCol$klDivergence, byRow$klDivergence, tolerance = 1e-10)
})

## ---- automatic k and method selection ---------------------------------------

test_that("k defaults to a valid value chosen by the eigengap heuristic", {
  mc <- randomChain(8, seed = 9)
  result <- aggregateStates(mc)
  expect_true(result$k >= 2L && result$k <= 7L)
  expect_length(result$partition, result$k)
})

test_that("method = 'adaptive' resolves to a concrete method and follows the documented rule", {
  small <- randomChain(10, seed = 2)
  resultSmall <- aggregateStates(small, k = 3, method = "adaptive")
  expect_equal(resultSmall$method, "spectral-top-down") # n < 30

  set.seed(11)
  n <- 32
  M <- matrix(runif(n * n), n, n)
  M <- M / rowSums(M)
  nm <- paste0("s", seq_len(n))
  dimnames(M) <- list(nm, nm)
  big <- new("markovchain", states = nm, transitionMatrix = M)

  resultBigSmallK <- aggregateStates(big, k = 4, method = "adaptive") # 4/32 <= 0.3
  expect_equal(resultBigSmallK$method, "spectral-bottom-up")

  resultBigLargeK <- aggregateStates(big, k = 28, method = "adaptive") # 28/32 > 0.3
  expect_equal(resultBigLargeK$method, "spectral-top-down")
})

## ---- validation and edge cases -----------------------------------------------

test_that("aggregateStates validates k", {
  mc <- randomChain(5, seed = 4)
  expect_error(aggregateStates(mc, k = 1), "between 2 and")
  expect_error(aggregateStates(mc, k = 5), "between 2 and")
  expect_error(aggregateStates(mc, k = 2.5), "single integer")
  expect_error(aggregateStates(mc, k = c(2, 3)), "single integer")
})

test_that("aggregateStates rejects reducible and periodic chains", {
  reducible <- new("markovchain", states = c("a", "b"),
                   transitionMatrix = matrix(c(1, 0, 0.5, 0.5),
                                              byrow = TRUE, nrow = 2,
                                              dimnames = list(c("a", "b"), c("a", "b"))))
  expect_error(aggregateStates(reducible, k = 1), "irreducible")

  cyc <- new("markovchain", states = c("a", "b", "c"),
             transitionMatrix = matrix(c(0, 1, 0,
                                         0, 0, 1,
                                         1, 0, 0), byrow = TRUE, nrow = 3,
                                       dimnames = list(c("a", "b", "c"), c("a", "b", "c"))))
  expect_error(aggregateStates(cyc, k = 2), "aperiodic")
})

test_that("aggregateStates requires at least 3 states", {
  two <- randomChain(2, seed = 6)
  expect_error(aggregateStates(two, k = 2), "at least 3 states")
})

## ---- defensive / low-level coverage -------------------------------------

test_that("aggregateStates rejects a corrupted (non-finite) transition matrix", {
  corrupted <- mcLump
  Pbad <- as.matrix(mcLump@transitionMatrix)
  Pbad["a", "b"] <- Inf
  corrupted@transitionMatrix <- Pbad
  expect_error(aggregateStates(corrupted, k = 2), "square and finite")
})
