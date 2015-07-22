library(markovchain)

statesNames=c("a","b","c")
markovB<-new("markovchain", states=statesNames, transitionMatrix=
               matrix(c(0.2,0.5,0.3,
                        0,1,0,
                        0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(statesNames,statesNames)
               ))

sequence<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
mcFit<-markovchainFit(data=sequence,byrow=FALSE)
# verifyMarkovProperty(mcFit)
verifyMarkovProperty(sequence)

# verifyMarkovProperty(markovB)
assessOrder(markovB)
assessStationarity(object)
divergenceTest(object)

test_that("States are those that should be", {
  # expect_equal(verifyMarkovProperty(markovB), TRUE)
})
