# TODO: Add comment
# 
# Author: Giorgio Spedicato
###############################################################################


#creates some markovchain objects
statesNames=c("a","b")
mcA<-new("markovchain", transitionMatrix=matrix(c(0.7,0.3,0.1,0.9),byrow=TRUE, nrow=2, 
				dimnames=list(statesNames,statesNames)
		))
mcB<-new("markovchain", states=c("a","b","c"), transitionMatrix=
				matrix(c(0.2,0.5,0.3,
								0,1,0,
								0.1,0.8,0.1),nrow=3, byrow=TRUE))
mcC<-new("markovchain", states=c("a","b","c","d"), 
		transitionMatrix=matrix(c(0.25,0.75,0,0,0.4,
						0.6,0,0,0,0,0.1,0.9,0,0,0.7,0.3), nrow=4, byrow=TRUE)
)
mcD<-new("markovchain", transitionMatrix=matrix(c(0,1,0,1), nrow=2,byrow=TRUE))

#apply some methods

testConversion<-as(mcC, "data.frame")
markovD<-t(mcC)
steadyStates(mcC)
steadyStates(mcC)

#perform some fitting

sequence<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
mcFit<-markovchainFit(data=sequence,byrow=FALSE)


