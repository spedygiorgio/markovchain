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
