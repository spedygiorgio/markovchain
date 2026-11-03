# Tests of the session cache of eigenvalues (R/eigenCache.R).

.randomChain <- function(n, seed, byrow = TRUE, name = "s") {
  set.seed(seed)
  P <- matrix(stats::rexp(n * n)^2, n)
  P <- P / rowSums(P)
  if (!byrow) P <- t(P)
  new("markovchain", transitionMatrix = P, states = paste0(name, seq_len(n)),
      byrow = byrow)
}

.cacheEntries <- function() length(markovchain:::.eigenCache$entries)

test_that("spectral functions give identical results with and without the cache", {
  old <- options(markovchain.cache = TRUE)
  on.exit(options(old))
  chains <- list(.randomChain(6, 1), .randomChain(9, 2, byrow = FALSE),
                 .randomChain(25, 3))
  # a periodic chain (complex eigenvalues on the unit circle)
  chains[[4]] <- new("markovchain", states = c("a", "b", "c"),
                     transitionMatrix = matrix(c(0, 1, 0, 0, 0, 1, 1, 0, 0), 3, byrow = TRUE))
  run <- function(mc) list(slem(mc), spectralGap(mc), relaxationTime(mc),
                           impliedTimescales(mc), topologicalEntropy(mc))
  for (mc in chains) {
    markovchain:::.clearEigenCache()
    first <- run(mc)    # fills the cache
    second <- run(mc)   # reads it
    options(markovchain.cache = FALSE)
    uncached <- run(mc)
    options(markovchain.cache = TRUE)
    expect_identical(first, uncached)
    expect_identical(second, uncached)
  }
})

test_that("the cached values are exactly those of eigen()", {
  markovchain:::.clearEigenCache()
  M <- .randomChain(12, 4)@transitionMatrix
  expect_identical(markovchain:::.eigenvaluesCached(M),
                   eigen(unname(M), only.values = TRUE)$values)
  expect_identical(markovchain:::.eigenvaluesCached(M),
                   eigen(unname(M), only.values = TRUE)$values)
})

test_that("a repeated matrix is a hit, a different one a miss", {
  old <- options(markovchain.cache = TRUE, markovchain.cacheSize = NULL)
  on.exit(options(old))
  markovchain:::.clearEigenCache()
  mc <- .randomChain(8, 5)
  slem(mc)
  expect_equal(.cacheEntries(), 1L)
  spectralGap(mc); relaxationTime(mc); impliedTimescales(mc)
  expect_equal(.cacheEntries(), 1L)

  # the same matrix with other state names: still a hit
  renamed <- .randomChain(8, 5, name = "other")
  slem(renamed)
  expect_equal(.cacheEntries(), 1L)

  # changing one entry by one ulp makes a different matrix: a new entry, and the value
  # must be the one of the new matrix, not the stored one
  M <- mc@transitionMatrix
  M2 <- M
  M2[1, 1] <- M2[1, 1] * (1 + .Machine$double.eps)
  v2 <- markovchain:::.eigenvaluesCached(M2)
  expect_equal(.cacheEntries(), 2L)
  expect_identical(v2, eigen(unname(M2), only.values = TRUE)$values)
})

test_that("least recently used entries are dropped to respect the budget", {
  old <- options(markovchain.cache = TRUE)
  on.exit(options(old))
  markovchain:::.clearEigenCache()
  Ms <- lapply(1:3, function(s) unname(.randomChain(20, 10 + s)@transitionMatrix))
  one <- as.numeric(utils::object.size(Ms[[1]])) +
    as.numeric(utils::object.size(eigen(Ms[[1]], only.values = TRUE)$values))
  options(markovchain.cacheSize = 2.5 * one)   # room for two entries

  markovchain:::.eigenvaluesCached(Ms[[1]])
  markovchain:::.eigenvaluesCached(Ms[[2]])
  markovchain:::.eigenvaluesCached(Ms[[1]])   # 1 becomes the most recently used
  markovchain:::.eigenvaluesCached(Ms[[3]])   # 2 is dropped
  stored <- lapply(markovchain:::.eigenCache$entries, `[[`, "matrix")
  expect_length(stored, 2L)
  expect_identical(stored[[1]], Ms[[3]])
  expect_identical(stored[[2]], Ms[[1]])

  # a matrix larger than the budget is not stored at all
  options(markovchain.cacheSize = 10)
  markovchain:::.clearEigenCache()
  markovchain:::.eigenvaluesCached(Ms[[1]])
  expect_equal(.cacheEntries(), 0L)

  # an invalid budget falls back to the default
  options(markovchain.cacheSize = "a lot")
  markovchain:::.eigenvaluesCached(Ms[[1]])
  expect_equal(.cacheEntries(), 1L)
})

test_that("options(markovchain.cache = FALSE) switches the cache off", {
  old <- options(markovchain.cache = FALSE)
  on.exit(options(old))
  markovchain:::.clearEigenCache()
  mc <- .randomChain(7, 20)
  slem(mc); spectralGap(mc)
  expect_equal(.cacheEntries(), 0L)
})

test_that("the cache does not change aggregateStates(k = NULL) or the eigenvalue plot", {
  old <- options(markovchain.cache = TRUE)
  on.exit(options(old))
  mc <- .randomChain(10, 30)
  markovchain:::.clearEigenCache()
  a1 <- aggregateStates(mc)
  options(markovchain.cache = FALSE)
  a2 <- aggregateStates(mc)
  expect_identical(a1, a2)

  skip_if_not_installed("ggplot2")
  options(markovchain.cache = TRUE)
  p1 <- ggplot2::autoplot(mc, type = "eigenvalues")
  options(markovchain.cache = FALSE)
  p2 <- ggplot2::autoplot(mc, type = "eigenvalues")
  expect_identical(p1$layers[[4]]$data, p2$layers[[4]]$data)
})
