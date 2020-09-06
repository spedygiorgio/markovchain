context("Checking verifyEmpiricalToTheoretical")

example_results <- list(statistic = 6.551795, dof = 6, pvalue = 0.3642899)

test_that("Sequence data input is computed correctly", {
  sequence <- c(0,1,2,2,1,0,0,0,0,0,0,1,2,2,2,1,0,0,1,0,0,0,0,0,0,1,1,
  2,0,0,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,2,1,0,
  0,2,1,0,0,0,0,0,0,1,1,1,2,2,0,0,2,1,1,1,1,2,1,1,1,1,1,1,1,1,1,0,2,
  0,1,1,0,0,0,1,2,2,0,0,0,0,0,0,2,2,2,1,1,1,1,0,1,1,1,1,0,0,2,1,1,
  0,0,0,0,0,2,2,1,1,1,1,1,2,1,2,0,0,0,1,2,2,2,0,0,0,1,1)
   
  mc <- matrix(c(5/8,1/4,1/8,1/4,1/2,1/4,1/4,3/8,3/8),
               byrow = TRUE,
               nrow = 3)
  rownames(mc) <- colnames(mc) <- 0:2
  theoreticalMc <- as(mc, "markovchain")
   
  result <- verifyEmpiricalToTheoretical(data = sequence, object = theoreticalMc, verbose = FALSE)
  
  expect_equivalent(example_results, result)
})

test_that("Matrix data input is computed correctly", {
  matrix <- matrix(c(51, 11, 8,
                     12, 31, 9,
                     6, 11, 10),
                   byrow = TRUE,
                   nrow = 3)
  rownames(matrix) <- colnames(matrix) <- 0:2
  
  mc <- matrix(c(5/8,1/4,1/8,1/4,1/2,1/4,1/4,3/8,3/8),
               byrow = TRUE,
               nrow = 3)
  rownames(mc) <- colnames(mc) <- 0:2
  theoreticalMc <- as(mc, "markovchain")
  
  result <- verifyEmpiricalToTheoretical(data = matrix, object = theoreticalMc, verbose = FALSE)
  
  expect_equivalent(example_results, result)
})

test_that("Input data sequences can contain missing states", {
  mc <- matrix(c(
    1 / 10, 7 / 10, 1 / 10, 1 / 10,
    1 / 10, 1 / 10, 4 / 10, 4 / 10,
    1 / 10, 5 / 10, 1 / 10, 3 / 10,
    1 / 10, 5 / 10, 3 / 10, 1 / 10
  ),
  byrow = TRUE,
  nrow = 4
  )
  rownames(mc) <- c(1:4)
  colnames(mc) <- c(1:4)
  theoreticalMc <- as(mc, "markovchain")
  
  sequence <- c(1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 4, 1, 1, 1, 1, 1, 1, 1)
  
  result <- verifyEmpiricalToTheoretical(data = sequence, object = theoreticalMc, verbose = FALSE)
  
  expect_true(exists("result"))
})

#TODO: review this part @DanielEbbert
# test_that("Null hypothesis is rejected when 0 in the object matrix is not 0 in the data matrix", {
#   example <- matrix(c(
#     0.6105, 0.1665, 0.0393, 0.1837,
#     0.1374, 0.5647, 0.0637, 0.2342,
#     0.3010, 0.1142, 0.3218, 0.2630,
#     0.2595, 0.3109, 0.0000, 0.4296
#   ),
#   byrow = TRUE,
#   nrow = 4
#   )
#   rownames(example) <- c(1:4)
#   colnames(example) <- c(1:4)
#   
#   mc <- matrix(c(
#     0.00, 1.00, 0.00, 0.00,
#     0.00, 0.00, 0.50, 0.50,
#     0.00, 0.75, 0.00, 0.25,
#     0.00, 0.75, 0.25, 0.00
#   ),
#   byrow = TRUE,
#   nrow = 4
#   )
#   rownames(mc) <- c(1:4)
#   colnames(mc) <- c(1:4)
#   theoreticalMc <- as(mc, "markovchain")
#   
#   result <- verifyEmpiricalToTheoretical(data = example, object = theoreticalMc, verbose = FALSE)
#   
#   expect_equivalent(result$pvalue, 0)
# })