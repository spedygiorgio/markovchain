library(markovchain)
library(MultinomialCI)

seq<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
mcfit<-markovchainFit(data=seq,byrow=TRUE)
# print(mcfit)
seqmat<-createSequenceMatrix(seq)
seqmat
mCI <- .multinomialCIRcpp(mcfit$estimate@transitionMatrix, seqmat, 0.95)
# print(mCI)

####end of creating multinomialCI
context("Multinomial confidence interval")

test_that("multinomial CI statisfay", {
#   expect_equal(mCI$lowerEndpointMatrix, matrix(c(0.2222222,0.3333333,
#                                                  0.5714286,0.1428571),nrow=2, byrow=TRUE, dimnames=list(c("a","b"),
#                                                                                                 c("a","b"))
#   ))
#   expect_equal(mCI$upperEndpointMatrix, matrix(c(0.8070205,0.9181316,
#                                                  1,0.6806468),nrow=2, byrow=TRUE, dimnames=list(c("a","b"),
#                                                                                                 c("a","b"))
#   ))
  expect_equal(mCI$upperEndpointMatrix[2,1],1)
})

# Multinomial distribution with 3 classes, from which 79 samples
# were drawn: 23 of them belong to the first class, 12 to the
# second class and 44 to the third class. Punctual estimations
# of the probabilities from this sample would be 23/79, 12/79
# and 44/79 but we want to build 95% simultaneous confidence intervals
# for the true probabilities
# m = multinomialCI(c(23,12,44), 0.05)
# print(paste("First class: [", m[1,1], m[1,2], "]"))
# print(paste("Second class: [", m[2,1], m[2,2], "]"))
# print(paste("Third class: [", m[3,1], m[3,2], "]"))

# seq<-c(4, 5)
# m = multinomialCI(seq, 0.05)
# m
