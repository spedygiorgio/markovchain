library(markovchain)

statesNames=c("a","b","c")
markovB<-new("markovchain", states=statesNames, transitionMatrix=
               matrix(c(0.2,0.5,0.3,
                        0,1,0,
                        0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(statesNames,statesNames)
               ))

sequence<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
mcFit<-markovchainFit(data=sequence,byrow=FALSE)
# verifyMarkovProperty(sequence)
# assessOrder(sequence)
assessStationarity(sequence)
divergenceTest(seqeuence)

data(blanden)
myMc<-as(blanden,"markovchain")
sequenza<-rmarkovchain(n = 100,myMc)
sequenza
res<-verifyMarkovProperty(sequenza)
res<-assessOrder(sequenza)
# print(res)

test_that("States are those that should be", {
  # expect_equal(verifyMarkovProperty(sequenza), TRUE)
})
