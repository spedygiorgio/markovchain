library(markovchain)

P <- matlab::zeros(10)
P[1, c(1, 3)] <- 1/2;
P[2, 2] <- 1/3; P[2,7] <- 2/3;
P[3, 1] <- 1;
P[4, 5] <- 1;
P[5, c(4, 5, 9)] <- 1/3;
P[6, 6] <- 1;
P[7, 7] <- 1/4; P[7,9] <- 3/4;
P[8, c(3, 4, 8, 10)] <- 1/4;
P[9, 2] <- 1;
P[10, c(2, 5, 10)] <- 1/3;
rownames(P) <- letters[1:10]
colnames(P) <- letters[1:10]
probMc <- new("markovchain", transitionMatrix = P,
             name = "Probability MC")

context("Recurrent Classes")

test_that("Checking recurrent classes for a known matrix", {
  expect_equal(.recurrentClassesRcpp(probMc), list(c("a", "c")
                                            , c("b", "g", "i")
                                            , c("f")))
})

