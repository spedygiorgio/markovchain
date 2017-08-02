context("Checking that fitting works")

#load & prepare data
data(rain)
data(holson)

myHolson<-as.matrix(holson[,-1]); rownames(myHolson)<-holson$id


test_that("Check createSequenceMatrix", {
  expect_equal(createSequenceMatrix(rain$rain), checksAlofiRawTransitions)
})

#data preparation
ciao<-c("a","a","b","b","a",NA,"b","a","b","a","a")

test_that("Check markovchainFit & listFit", {
  expect_equal(markovchainFit(ciao), simpleMcCiaoFit)
  expect_equal(markovchainListFit(data=myHolson), checkmarkovchainFitList)
})


#### tests for noofVisitsDist function

transMatr<-matrix(c(0.4,0.6,.3,.7),nrow=2,byrow=TRUE)
simpleMc<-new("markovchain", states=c("a","b"),
              transitionMatrix=transMatr, 
              name="simpleMc")  

answer <- c(0.348148, 0.651852)
names(answer) <- c("a","b")
test_that("Check noofVisitsDist works", {
  expect_equal(noofVisitsDist(simpleMc,5,"a"),answer)
})

