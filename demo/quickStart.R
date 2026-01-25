# TODO: Add comment
# 
# Author: Giorgio Spedicato
###############################################################################
require(markovchain)

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


#canonic form

#from https://math.dartmouth.edu/archive/m20x06/public_html/Lecture14.pdf

P<-markovchain::zeros(5)
P[1,1]<-P[5,5]<-1
P[2,1]<-P[2,3]<-0.5
P[3,2]<-P[3,4]<-0.5
P[4,3]<-P[4,5]<-0.5
mcP<-as(P,"markovchain")
mcPCan<-canonicForm(mcP)


# coercing markov chains to sparse matrix forth and back
require(Matrix)
ciauz<-c(0,.5,.5,1,0,0,.2,0,.8)
matrix(ciauz, nrow=3, byrow=TRUE)
sparse<-as(matrix(ciauz, nrow=3, byrow=TRUE),"sparseMatrix")
mc<-as(sparse,"markovchain")
sparse2<-as(mc,"sparseMatrix")
