context("randomMarkovChain(), dirichletChain() and identityChain()")

test_that("randomMarkovChain returns a valid chain, reproducible with seed", {
  a <- randomMarkovChain(5, seed = 3)
  expect_s4_class(a, "markovchain")
  expect_equal(unname(rowSums(a@transitionMatrix)), rep(1, 5), tolerance = 1e-12)
  expect_identical(a@states, as.character(1:5))
  expect_identical(randomMarkovChain(5, seed = 3), a)
  expect_false(identical(randomMarkovChain(5, seed = 4)@transitionMatrix,
                         a@transitionMatrix))
})

test_that("seed does not touch the caller's random stream", {
  set.seed(10); before <- runif(1)
  set.seed(10); randomMarkovChain(4, seed = 99); dirichletChain(4, 2, seed = 5)
  expect_identical(runif(1), before)
  # without seed, set.seed() beforehand controls the result
  set.seed(1); x <- randomMarkovChain(3)
  set.seed(1); y <- randomMarkovChain(3)
  expect_identical(x, y)
})

test_that("zeros gives exactly that many zero probabilities", {
  for (z in c(0, 3, 12)) {
    P <- randomMarkovChain(4, zeros = z, seed = z)@transitionMatrix
    expect_equal(sum(P == 0), z)
    expect_true(all(rowSums(P > 0) >= 1))
  }
  expect_error(randomMarkovChain(4, zeros = 13), "cannot exceed 12")
})

test_that("mask fixes probabilities and the free entries share the rest", {
  m <- matrix(NA, 3, 3)
  m[1, 2] <- 0.3
  m[2, ] <- c(0.5, 0.5, NA)      # full row: the NA becomes 0
  P <- randomMarkovChain(states = c("a", "b", "c"), mask = m, seed = 1)@transitionMatrix
  expect_equal(P["a", "b"], 0.3)
  expect_equal(unname(P["b", ]), c(0.5, 0.5, 0))
  expect_equal(unname(rowSums(P)), rep(1, 3), tolerance = 1e-12)
  # byrow = FALSE reads the mask by columns
  Q <- randomMarkovChain(3, mask = t(m), byrow = FALSE, seed = 1)
  expect_false(Q@byrow)
  expect_equal(unname(Q@transitionMatrix[, 2]), c(0.5, 0.5, 0))
  expect_equal(unname(colSums(Q@transitionMatrix)), rep(1, 3), tolerance = 1e-12)
})

test_that("randomMarkovChain validates its arguments", {
  expect_error(randomMarkovChain(0), "at least 1")
  expect_error(randomMarkovChain(3, states = c("a", "a", "b")), "distinct")
  expect_error(randomMarkovChain(3, zeros = -1), "zeros")
  bad <- matrix(NA, 2, 2); bad[1, ] <- c(0.7, 0.7)
  expect_error(randomMarkovChain(2, mask = bad), "more than one")
  bad <- matrix(NA, 2, 2); bad[1, ] <- c(0.2, 0.2)
  expect_error(randomMarkovChain(2, mask = bad), "no entry is left free")
  expect_error(randomMarkovChain(2, mask = matrix(NA, 3, 3)), "2 x 2")
  expect_error(randomMarkovChain(2, seed = 1.5), "seed")
})

test_that("zeros are spread uniformly over the free entries", {
  # Every free entry not reserved to be positive is zero with probability
  # zeros / (#free - #rows), times the chance of not being the reserved one.
  set.seed(2024)
  N <- 3000
  counts <- matrix(0, 3, 3)
  for (s in seq_len(N)) counts <- counts + (randomMarkovChain(3, zeros = 3)@transitionMatrix == 0)
  expected <- (2 / 3) * (3 / 6)
  expect_equal(as.numeric(counts / N), rep(expected, 9), tolerance = 0.04)
})

test_that("dirichletChain follows the stick-breaking construction", {
  d <- dirichletChain(6, diffusion = 2, seed = 1)
  expect_equal(unname(rowSums(d@transitionMatrix)), rep(1, 6), tolerance = 1e-12)
  expect_identical(dirichletChain(6, diffusion = 2, seed = 1), d)
  # E[w_1] before normalisation is 1 / (1 + alpha): the first column carries
  # the largest average probability, the last the smallest
  set.seed(5)
  M <- Reduce(`+`, lapply(1:400, function(i) dirichletChain(5, 1)@transitionMatrix)) / 400
  cm <- colMeans(M)
  expect_true(cm[1] > cm[5])
  # shiftConcentration reverses the columns
  set.seed(5)
  Ms <- Reduce(`+`, lapply(1:400, function(i)
    dirichletChain(5, 1, shiftConcentration = TRUE)@transitionMatrix)) / 400
  expect_equal(unname(colMeans(Ms)), unname(rev(cm)), tolerance = 0.05)
  # diagonalBias raises the diagonal
  set.seed(6)
  B <- Reduce(`+`, lapply(1:400, function(i)
    dirichletChain(5, 3, diagonalBias = 5)@transitionMatrix)) / 400
  expect_true(mean(diag(B)) > mean(B[row(B) != col(B)]))
  expect_error(dirichletChain(1, 1), "at least 2")
  expect_error(dirichletChain(3, 0), "diffusion")
  expect_error(dirichletChain(3, 1, diagonalBias = -1), "diagonalBias")
})

test_that("identityChain is the identity", {
  i3 <- identityChain(3)
  expect_equal(unname(i3@transitionMatrix), diag(3))
  expect_identical(absorbingStates(i3), as.character(1:3))
  i2 <- identityChain(states = c("x", "y"))
  expect_identical(i2@states, c("x", "y"))
  mc <- randomMarkovChain(3, seed = 2)
  expect_equal((identityChain(3) * mc)@transitionMatrix, mc@transitionMatrix)
  expect_error(identityChain(), "Provide n or states")
})
