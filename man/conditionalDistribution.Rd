\name{conditionalDistribution}
\alias{conditionalDistribution}
\title{\code{conditionalDistribution} of a Markov Chain}
\description{
It extracts the conditional distribution of the subsequent state, given current state.
}
\usage{
conditionalDistribution(object,state)
}
\arguments{
  \item{object}{A \code{markovchain} object.}
  \item{state}{Subsequent state.}
}
\value{A named probability vector}
\references{A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall 2010}
\author{
Giorgio Spedicato
}
\seealso{
\code{\linkS4class{markovchain}}
}
\examples{
#define a markov chain
statesNames=c("a","b","c")
markovB<-new("markovchain", states=statesNames, transitionMatrix=matrix(c(0.2,0.5,0.3,
0,1,0,0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(statesNames,statesNames)))
conditionalDistribution(markovB,"b")
}


