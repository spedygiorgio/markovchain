#library(markovchain)

sequence<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
mcFit<-markovchainFit(data=sequence,byrow=FALSE)
# verifyMarkovProperty(sequence)
# assessOrder(sequence)
# assessStationarity(sequence, 1)
# divergenceTest(sequence, mcFit$estimate@transitionMatrix)

data(blanden)
myMc<-as(blanden,"markovchain")
# print(myMc)
sequenza<-rmarkovchain(n = 100,myMc)
sequenza
res<-verifyMarkovProperty(sequenza)
res<-assessOrder(sequenza)
res<-assessStationarity(sequenza, 10)
res<-divergenceTest(sequenza, myMc)
# print(res)

test_that("States are those that should be", {
  # expect_equal(verifyMarkovProperty(sequenza), TRUE)
})
