context("Lumpability of Markov Chains")

test_that("Helper function strictly validates partitions", {
  mc <- new("markovchain", states = c("A", "B", "C"), 
            transitionMatrix = matrix(1/3, 3, 3, dimnames=list(c("A","B","C"), c("A","B","C"))))
  
  # Missing a state (C is missing)
  bad_part_1 <- list(Group1 = c("A", "B"))
  expect_error(lump(mc, bad_part_1, force = TRUE), "no omissions")
  
  # Duplicated state (B is in both)
  bad_part_2 <- list(Group1 = c("A", "B"), Group2 = c("B", "C"))
  expect_error(lump(mc, bad_part_2, force = TRUE), "no duplicates")
  
  # Typo in state name
  bad_part_3 <- list(Group1 = c("A", "B"), Group2 = c("Z"))
  expect_error(lump(mc, bad_part_3, force = TRUE), "do not exist in the Markov chain")
})


test_that("is.lumpable and lump work correctly on 'Land of Oz' textbook example", {
  # The classic Grinstead & Snell (1997) "Land of Oz" weather model
  states_oz <- c("rainy", "nice", "snowy")
  P_oz <- matrix(c(0.5,  0.25, 0.25,
                   0.5,  0.0,  0.5,
                   0.25, 0.25, 0.5), 
                 byrow = TRUE, nrow = 3, dimnames = list(states_oz, states_oz))
  
  mcOz <- new("markovchain", states = states_oz, transitionMatrix = P_oz, name = "Land of Oz")
  
  # We partition the weather into "Bad" (Rainy/Snowy) and "Nice"
  partition_oz <- list(
    Bad_Weather = c("rainy", "snowy"),
    Nice_Weather = c("nice")
  )
  
  # 1. The chain is exactly lumpable under this partition
  expect_true(is.lumpable(mcOz, partition_oz))
  
  # 2. Verify the output lumped Markov chain
  lumped_oz <- lump(mcOz, partition_oz)
  expect_s4_class(lumped_oz, "markovchain")
  expect_equal(states(lumped_oz), c("Bad_Weather", "Nice_Weather"))
  
  # In Land of Oz, probability of going from Bad to Nice is ALWAYS 0.25
  expect_equal(lumped_oz@transitionMatrix["Bad_Weather", "Nice_Weather"], 0.25)
  # Probability of staying in Bad Weather is 0.75
  expect_equal(lumped_oz@transitionMatrix["Bad_Weather", "Bad_Weather"], 0.75)
})


test_that("is.lumpable and lump work correctly on empirical Jukes-Cantor DNA model", {
  # Jukes-Cantor (1969) model for DNA sequence evolution
  # Bases mutate to any other base with uniform probability alpha.
  alpha <- 0.05
  P_jc <- matrix(alpha, nrow = 4, ncol = 4)
  diag(P_jc) <- 1 - 3 * alpha
  states_dna <- c("A", "C", "G", "T")
  rownames(P_jc) <- colnames(P_jc) <- states_dna
  
  mcJC <- new("markovchain", states = states_dna, transitionMatrix = P_jc)
  
  # Chemical partition: Purines (A, G) and Pyrimidines (C, T)
  partition_dna <- list(
    Purines = c("A", "G"),
    Pyrimidines = c("C", "T")
  )
  
  # The uniform symmetry makes this exactly lumpable
  expect_true(is.lumpable(mcJC, partition_dna))
  
  lumped_dna <- lump(mcJC, partition_dna)
  
  # Cumulative probability of mutating from a Purine to ANY Pyrimidine is 2 * alpha
  expect_equal(lumped_dna@transitionMatrix["Purines", "Pyrimidines"], 2 * alpha)
})


test_that("lump function correctly stops on non-lumpable matrices without force = TRUE", {
  # A completely asymmetrical, non-lumpable matrix
  P_bad <- matrix(c(0.5, 0.3, 0.2,
                    0.7, 0.2, 0.1,
                    0.1, 0.4, 0.5), 
                  byrow = TRUE, nrow = 3, dimnames = list(c("A1","A2","B"), c("A1","A2","B")))
  mcBad <- new("markovchain", states = c("A1", "A2", "B"), transitionMatrix = P_bad)
  
  partition_bad <- list(Group_A = c("A1", "A2"), Group_B = c("B"))
  
  # Fails the exact lumpability test
  expect_false(is.lumpable(mcBad, partition_bad))
  
  # Should throw an error preventing mathematically invalid aggregation
  expect_error(lump(mcBad, partition_bad), "not exactly lumpable")
  
  # Should succeed if the user consciously forces the aggregation
  expect_s4_class(lump(mcBad, partition_bad, force = TRUE), "markovchain")
})


test_that("autoLump performs approximate Spectral Clustering on S&P real data", {
  # Standard & Poor's empirical rating transition matrix (Real world data)
  rc <- c("AAA", "AA", "A", "BBB", "BB", "B", "CCC", "D")
  creditMatrix <- matrix(
    c(90.81, 8.33, 0.68, 0.06, 0.08, 0.02, 0.01, 0.01,
      0.70, 90.65, 7.79, 0.64, 0.06, 0.13, 0.02, 0.01,
      0.09, 2.27, 91.05, 5.52, 0.74, 0.26, 0.01, 0.06,
      0.02, 0.33, 5.95, 85.93, 5.30, 1.17, 1.12, 0.18,
      0.03, 0.14, 0.67, 7.73, 80.53, 8.84, 1.00, 1.06,
      0.01, 0.11, 0.24, 0.43, 6.48, 83.46, 4.07, 5.20,
      0.21, 0, 0.22, 1.30, 2.38, 11.24, 64.86, 19.79,
      0, 0, 0, 0, 0, 0, 0, 100
    )/100, 8, 8, dimnames = list(rc, rc), byrow = TRUE)
  
  creditMc <- new("markovchain", states = rc, transitionMatrix = creditMatrix)
  
  # Automatically discover 3 latent macro-states
  k <- 3
  result <- autoLump(creditMc, k)
  
  expect_type(result, "list")
  expect_true(!is.null(result$partition))
  expect_s4_class(result$lumped_chain, "markovchain")
  
  # The resulting lumped chain must have exactly 3 dimensions
  expect_equal(length(states(result$lumped_chain)), k)
  expect_equal(dim(result$lumped_chain@transitionMatrix), c(k, k))
})