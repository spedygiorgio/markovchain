# Period examples from http://www.math.wisc.edu/~anderson/605F11/Notes/StochBioChapter3.pdf

#library(markovchain)

mcPeriodic<-new("markovchain", states=c("0","1","2"), transitionMatrix=
                 matrix(c(0,1,0,
                          0,0,1,
                          1,0,0),nrow=3, byrow=TRUE, dimnames=list(c("0","1","2"),
                                                                   c("0","1","2"))
                 ))

mcAperiodic<-new("markovchain", states=c("0","1","2","3","4"), transitionMatrix=
                 matrix(c(1/2,1/2,0,0,0,
                          1/2,0,1/2,0,0,
                          0,1/2,0,1/2,0,
                          0,0,1/2,0,1/2,
                          0,0,0,1,0),nrow=5, byrow=TRUE, dimnames=list(c("0","1","2","3","4"),
                                                                   c("0","1","2","3","4"))
                 ))

####end creating DTMC
context("Basic DTMC proprieties")

test_that("States are those that should be", {
  expect_equal(period(mcPeriodic),3)
  expect_equal(period(mcAperiodic),1)
})