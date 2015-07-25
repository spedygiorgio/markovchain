\name{verifyMarkovProperty}
\alias{verifyMarkovProperty}
\alias{assessOrder}
\alias{assessStationarity}
\alias{divergenceTest}
\title{
Various functions to perform statistical inference of DTMC
}
\description{
These functions verify the Markov property, assess the order and stationarity of the \code{markovchain} objects.
}
\usage{
verifyMarkovProperty(mc)
assessOrder(mc)
assessStationarity(object)
divergenceTest(object)
}
\arguments{
  \item{mc}{A markov chain sequence.}
  \item{object}{A \code{markovchain} object.}
}
\value{Verification result or the order}
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
sequence<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
verifyMarkovProperty(sequence)
assessOrder(markovB)
assessStationarity(object)
divergenceTest(object)
}

