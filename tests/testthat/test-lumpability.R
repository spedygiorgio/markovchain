context("Lumpability of Markov Chains")

test_that("helper function strictly validates partitions", {
  mc <- new("markovchain", states = c("A", "B", "C"),
            transitionMatrix = matrix(1 / 3, 3, 3,
                                      dimnames = list(c("A", "B", "C"), c("A", "B", "C"))))

  expect_error(lump(mc, list(Group1 = c("A", "B")), force = TRUE), "no duplicates, no omissions")
  expect_error(lump(mc, list(Group1 = c("A", "B"), Group2 = c("B", "C")), force = TRUE), "no duplicates, no omissions")
  expect_error(lump(mc, list(Group1 = c("A", "B"), Group2 = c("Z")), force = TRUE), "do not exist")
})

test_that("is.lumpable and lump work on the Land of Oz example", {
  states_oz <- c("rainy", "nice", "snowy")
  P_oz <- matrix(c(0.5,  0.25, 0.25,
                   0.5,  0.0,  0.5,
                   0.25, 0.25, 0.5),
                 byrow = TRUE, nrow = 3, dimnames = list(states_oz, states_oz))
  mcOz <- new("markovchain", states = states_oz, transitionMatrix = P_oz, name = "Land of Oz")
  partition_oz <- list(Bad_Weather = c("rainy", "snowy"), Nice_Weather = c("nice"))

  expect_true(is.lumpable(mcOz, partition_oz))
  lumped_oz <- lump(mcOz, partition_oz)
  expect_s4_class(lumped_oz, "markovchain")
  expect_equal(states(lumped_oz), c("Bad_Weather", "Nice_Weather"))
  expect_equal(lumped_oz@transitionMatrix["Bad_Weather", "Nice_Weather"], 0.25)
  expect_equal(lumped_oz@transitionMatrix["Bad_Weather", "Bad_Weather"], 0.75)
})

test_that("is.lumpable and lump work on the Jukes-Cantor DNA model", {
  alpha <- 0.05
  states_dna <- c("A", "C", "G", "T")
  P_jc <- matrix(alpha, nrow = 4, ncol = 4, dimnames = list(states_dna, states_dna))
  diag(P_jc) <- 1 - 3 * alpha
  mcJC <- new("markovchain", states = states_dna, transitionMatrix = P_jc)
  partition_dna <- list(Purines = c("A", "G"), Pyrimidines = c("C", "T"))

  expect_true(is.lumpable(mcJC, partition_dna))
  lumped_dna <- lump(mcJC, partition_dna)
  expect_equal(lumped_dna@transitionMatrix["Purines", "Pyrimidines"], 2 * alpha)
})

test_that("lump stops on non-lumpable matrices unless force = TRUE", {
  P_bad <- matrix(c(0.5, 0.3, 0.2,
                    0.7, 0.2, 0.1,
                    0.1, 0.4, 0.5),
                  byrow = TRUE, nrow = 3, dimnames = list(c("A1", "A2", "B"), c("A1", "A2", "B")))
  mcBad <- new("markovchain", states = c("A1", "A2", "B"), transitionMatrix = P_bad)
  partition_bad <- list(Group_A = c("A1", "A2"), Group_B = c("B"))

  expect_false(is.lumpable(mcBad, partition_bad))
  expect_error(lump(mcBad, partition_bad), "not exactly lumpable")
  expect_s4_class(lump(mcBad, partition_bad, force = TRUE), "markovchain")
})

test_that("autoLump performs approximate spectral clustering on credit-rating data", {
  rc <- c("AAA", "AA", "A", "BBB", "BB", "B", "CCC", "D")
  creditMatrix <- matrix(
    c(90.81, 8.33, 0.68, 0.06, 0.08, 0.02, 0.01, 0.01,
      0.70, 90.65, 7.79, 0.64, 0.06, 0.13, 0.02, 0.01,
      0.09, 2.27, 91.05, 5.52, 0.74, 0.26, 0.01, 0.06,
      0.02, 0.33, 5.95, 85.93, 5.30, 1.17, 1.12, 0.18,
      0.03, 0.14, 0.67, 7.73, 80.53, 8.84, 1.00, 1.06,
      0.01, 0.11, 0.24, 0.43, 6.48, 83.46, 4.07, 5.20,
      0.21, 0, 0.22, 1.30, 2.38, 11.24, 64.86, 19.79,
      0, 0, 0, 0, 0, 0, 0, 100) / 100,
    8, 8, dimnames = list(rc, rc), byrow = TRUE)
  creditMc <- new("markovchain", states = rc, transitionMatrix = creditMatrix)

  k <- 3
  result <- autoLump(creditMc, k)
  expect_type(result, "list")
  expect_true(!is.null(result$partition))
  expect_s4_class(result$lumped_chain, "markovchain")
  expect_equal(length(states(result$lumped_chain)), k)
  expect_equal(dim(result$lumped_chain@transitionMatrix), c(k, k))
})
