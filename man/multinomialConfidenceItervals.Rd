\name{multinomialConfidenceItervals}
\alias{multinomialConfidenceItervals}
\title{
Various function to perform structural analysis of DTMC
}
\description{
These functions return absorbing and transient states of the \code{markovchain} objects.
}
\usage{
multinomialConfidenceItervals(transitionMatrix, countsTransitionMatrix, confidencelevel)
}
\arguments{
  \item{transitionMatrix}{Transition matrix.}
  \item{countsTransitionMatrix}{Counts transition matrix.}
  \item{confidencelevel}{Confidence level.}
}
\value{A matrix}
\references{
Feres, Matlab listing for markov chain.
}
\author{
Tae Seung Kang, Giorgio Alfredo Spedicato
}
\seealso{
\code{\linkS4class{markovchain}}
}
\examples{
statesNames=c("a","b","c")
markovB<-new("markovchain", states=statesNames, transitionMatrix=
          matrix(c(0.2,0.5,0.3,
                   0,1,0,
                   0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(statesNames,statesNames)
                 ))
absorbingStates(markovB)
transientStates(markovB)
#periodicity analysis
E=matrix(0,nrow=4,ncol=4)
E[1,2]=1
#E[2,c(1, 3)]=1/2;
E[2,1]=1/3;E[2,3]=2/3
#E[3,c(2, 4)]=1/2;
E[3,2]=1/4;E[3,4]=3/4
E[4,3]=1
mcE <- new("markovchain", states = c("a", "b", "c", "d"), 
                 transitionMatrix = E, 
                 name = "E")

is.irreducible(mcE) #true
period(mcE) #2

myMatr<-matrix(c(0,0,1/2,1/4,1/4,0,0,
0,0,1/3,0,2/3,0,0,
0,0,0,0,0,1/3,2/3,
0,0,0,0,0,1/2,1/2,
0,0,0,0,0,3/4,1/4,
1/2,1/2,0,0,0,0,0,
1/4,3/4,0,0,0,0,0),byrow=TRUE,ncol=7)
myMc<-new("markovchain",transitionMatrix=myMatr)
period(myMc)
}

