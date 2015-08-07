\name{verifyMarkovProperty}
\alias{verifyMarkovProperty}
\alias{assessOrder}
\alias{assessStationarity}
\alias{divergenceTest}
\title{
Various functions to perform statistical inference of DTMC
}
\description{
These functions verify the Markov property, assess the order and stationarity of the Markov chain.
}
\usage{
verifyMarkovProperty(mc)
assessOrder(mc)
assessStationarity(mc)
divergenceTest(m1, m2, mc)
}
\arguments{
  \item{mc}{A markov chain sequence.}
  \item{m1}{A transition matrix for a markov chain sequence.}
  \item{m2}{A transition matrix for a markov chain sequence.}
}
\value{Verification result}
\references{
Monika, Anderson and Goodman.
}
\author{
Tae Seung Kang, Giorgio Alfredo Spedicato
}
\seealso{
\code{\linkS4class{markovchain}}
}
\examples{
sequence<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
mcFit<-markovchainFit(data=sequence,byrow=FALSE)
verifyMarkovProperty(sequence)
assessOrder(sequence)
assessStationarity(sequence)
divergenceTest(mcFit$estimate@transitionMatrix, mcFit$estimate@transitionMatrix, sequence)
}

