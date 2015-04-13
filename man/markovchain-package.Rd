\name{markovchain-package}
\alias{markovchain-package}
\alias{markovchain}
\docType{package}
\title{
Easy Handling Discrete Time Markov Chains
}
\description{
The package contains classes and method to create and manage (plot, print, export for example) discrete time Markov chains (DTMC). In addition it provide functions to perform statistical (fitting and drawing random variates) and probabilistic (analysis of DTMC proprieties) analysis.
}
\details{
\tabular{ll}{
Package: \tab markovchain\cr
Type: \tab Package\cr
Version: \tab 0.1.9\cr
Date: \tab 2015-04-15\cr
License: \tab GPL-2\cr
Depends: \tab R (>= 2.14), methods, expm, matlab, igraph, Matrix\cr
}
}
\author{
Giorgio Alfredo Spedicato
Maintainer: Giorgio Alfredo Spedicato <spedicato_giorgio@yahoo.it>

}
\references{
Discrete-Time Markov Models, Bremaud, Springer 1999
}
\keyword{ package }

\examples{
#create some markov chains
statesNames=c("a","b")
 mcA<-new("markovchain", transitionMatrix=matrix(c(0.7,0.3,0.1,0.9),byrow=TRUE,
 nrow=2, dimnames=list(statesNames,statesNames)))

statesNames=c("a","b","c")
mcB<-new("markovchain", states=statesNames, transitionMatrix=
          matrix(c(0.2,0.5,0.3,
                   0,1,0,
                   0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(statesNames,
				   statesNames)
                 ))

statesNames=c("a","b","c","d")
matrice<-matrix(c(0.25,0.75,0,0,0.4,0.6,0,0,0,0,0.1,0.9,0,0,0.7,0.3), 
nrow=4, byrow=TRUE)
mcC<-new("markovchain", states=statesNames, transitionMatrix=matrice)
mcD<-new("markovchain", transitionMatrix=matrix(c(0,1,0,1), nrow=2,byrow=TRUE))



#operations with S4 methods

mcA^2
steadyStates(mcB)
absorbingStates(mcB)
markovchainSequence(n=20, markovchain=mcC, include=TRUE)

}
