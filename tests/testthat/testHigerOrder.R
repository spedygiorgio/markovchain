library(markovchain)

sequence<-c("a", "a", "b", "b", "a", "c", "b", "a", "b", "c", "a", "b", "c", "a", "b", "c", "a", "b", "a", "b")
# mcFit<-fitHigherOrder(data=sequence,byrow=FALSE)
fitHigherOrder(sequence)
