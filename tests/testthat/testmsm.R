#library(markovchain)
library("msm")
context("Checking conversion of objects to msm")


# cav
cav <- msm::cav
Q <- rbind ( c(0, 0.25, 0, 0.25),
             c(0.166, 0, 0.166, 0.166),
             c(0, 0.25, 0, 0.25),
             c(0, 0, 0, 0) )
cavmsm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = Q, death = 4)
# cavmsm
# qmatrix.msm(cavmsm)
# qmatrix.msm(cavmsm, covariates = list(sex=1))
# pmatrix.msm(cavmsm)

# .msm2Mc(cavmsm)
msmMc <- as(cavmsm, "markovchain")
# print(msmMc)

test_that("Conversion of objects", 
          {
            expect_equal(class(msmMc)=="markovchain",TRUE)
          })