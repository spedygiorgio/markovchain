\name{transitionProbability}
\alias{transitionProbability}
\title{
Function to get the transition probabilities from initial to subsequent status.
}
\description{
This is a convenience function to get transition probabilities.
}
\usage{
transitionProbability(object,t0,t1)
}
\arguments{
  \item{object}{A \code{markovchain} object.}
  \item{t0}{Initial state.}
  \item{t1}{Subsequent state.}
}
\value{A matrix}
\references{
A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall 2010
}
\author{
Giorgio Spedicato
}
\seealso{
\code{\linkS4class{markovchain}}
}
\examples{
#define a markov chain
statesNames=c("a","b","c")
markovB<-new("markovchain", states=statesNames, transitionMatrix=
          matrix(c(0.2,0.5,0.3,
                   0,1,0,
                   0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(statesNames,statesNames)
                 ))
transitionProbability(markovB,"b", "c")
}

