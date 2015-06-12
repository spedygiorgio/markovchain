## ----setup, include=FALSE------------------------------------------------
library(knitr)
rm(list=ls())

## ----load, echo=TRUE, warning=FALSE--------------------------------------
library(markovchain) #load the package

## ----create, echo=TRUE, tidy=FALSE, message=FALSE,tidy=FALSE-------------
tmA <- matrix(c(0,0.5,0.5,.5,0,.5,.5,.5,0),nrow = 3,
              byrow = TRUE) #define the transition matrix
dtmcA <- new("markovchain",transitionMatrix=tmA, 
             states=c("a","b","c"), 
             name="MarkovChain A") #create the DTMC
dtmcA

## ----create2, echo=TRUE, tidy=FALSE--------------------------------------
dtmcA2<-as(tmA, "markovchain") #using coerce from matrix
states(dtmcA2) #note default names assigned to states

## ----plot, echo=TRUE-----------------------------------------------------
plot(dtmcA)

## ----proprieties, echo=TRUE, tidy=FALSE----------------------------------
dtmcA[2,3] #using [ method
transitionProbability(dtmcA, 
                      "b","c") #using specific S4 method
conditionalDistribution(dtmcA,"b")

## ----transitions, echo=TRUE, tidy=FALSE----------------------------------
initialState<-c(0,1,0)
steps<-4
finalState<-initialState*dtmcA^steps #using power operator
finalState

## ----steadystate, echo=TRUE, tidy=FALSE----------------------------------
steadyStates(dtmcA) #S4 method

## ----mathematicaMc, echo=TRUE, tidy=FALSE--------------------------------
E <- matrix(0, nrow = 4, ncol = 4)
E[1, 2] <- 1;E[2, 1] <- 1/3; E[2, 3] <- 2/3
E[3,2] <- 1/4; E[3, 4] <- 3/4; E[4, 3] <- 1
mcMathematica <- new("markovchain", states = c("a", "b", "c", "d"),
                     transitionMatrix = E,name = "Mathematica")

## ----summary, echo=TRUE, tidy=FALSE--------------------------------------
summary(mcMathematica)

## ----fitIntro, echo=TRUE, tidy=FALSE-------------------------------------
#using Alofi rainfall dataset
data(rain) 
mysequence<-rain$rain
createSequenceMatrix(mysequence)

## ----mlfFit, echo=TRUE, tidy=FALSE---------------------------------------
myFit<-markovchainFit(data=mysequence,confidencelevel = .9,method = "mle")
myFit

## ----mlfFit2, echo=TRUE, tidy=FALSE--------------------------------------
alofiMc<-myFit$estimate
alofiMc

