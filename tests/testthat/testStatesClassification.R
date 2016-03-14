context("Classification of states")

A <-  structure(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
              0.5, 0.3, 0.3, 0, 0, 0, 0, 0, 0.5, 0.7, 0.7, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0.3, 0.4, 0, 0, 0, 0, 0, 0.4, 0, 0.5, 0, 0, 0, 0, 0,
              0.6, 0.7, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 1), .Dim = c(8L, 8L), .Dimnames = list(
                c("1", "2", "3", "4", "5", "6", "7", "8"), c("1", "2", "3",
                                                             "4", "5", "6", "7", "8")))
mchain <- new("markovchain", transitionMatrix=A)

#summary(mchain)


test_that("States are those that should be", {
    expect_equal(recurrentClasses(mchain), list(c("3", "4"), c("8")))
    expect_equal(transientStates(mchain), c("1", "2", "5", "6", "7"))
    expect_equal(absorbingStates(mchain), "8")
})

#https://www.math.ucdavis.edu/~gravner/MAT135B/materials/ch13.pdf

mcMatr1<-matlab::zeros(3)
mcMatr1[1,]<-c(0.5,0.5,0)
mcMatr1[2,]<-c(0.5,0.25,0.25)
mcMatr1[3,]<-c(0,1/3,2/3)
mc1<-as(mcMatr1,"markovchain")

test_that("States are those that should be", {
expect_equal(is.irreducible(mc1),TRUE)
})

mcMatr2<-matrix(c(0, 0, 1/2, 1/2,1, 0 ,0, 0,0, 1, 0, 0,0, 1, 0, 0),ncol = 4,byrow=TRUE)
mc2<-as(mcMatr2,"markovchain")


test_that("States are those that should be", {
  expect_equal(recurrentClasses(mc2),list(c("s1","s2","s3","s4")))
})


mcMatr3<-matrix(c(
0,1,0,0,0,0,
0.4,0.6,0,0,0,0,
0.3,0,0.4,0.2,0.1,0,
0,0,0,0.3,0.7,0,
0,0,0,0.5,0,0.5,
0,0,0,0.3,0,0.7),nrow = 6,byrow=TRUE)

mc3<-as(mcMatr3,"markovchain")
recurrentClasses(mc3)
transientStates(mc3)
#canonicForm(mc3)


test_that("States are those that should be", {
  expect_equal(recurrentClasses(mc3),list(c("s1","s2"),c("s4","s5","s6") ))
  expect_equal(transientStates(mc3),"s3")
})

mcMatr4<-matlab::zeros(5)
mcMatr4[1:2,1:2]<-0.5*matlab::ones(2)
mcMatr4[5,1]<-1
mcMatr4[3,3]<-1
mcMatr4[4,3:4]<-0.5
mc4<-as(mcMatr4,"markovchain")


test_that("States are those that should be", {
  expect_equal(recurrentClasses(mc4),list(c("s1","s2"),c("s3")))
               expect_equal(absorbingStates(mc4),"s3")
               expect_equal(transientStates(mc4),c("s4","s5"))
  }
)
