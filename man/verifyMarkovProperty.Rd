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
verifyMarkovProperty(sequence, ...)
assessOrder(sequence)
assessStationarity(sequence, nblocks)
divergenceTest(sequence, hypothetic)
}
\arguments{
  \item{sequence}{An empirical sequence.}
  \item{...}{Parameters for chi-square test.}
  \item{hypothetic}{A transition matrix for a hypothetic markov chain sequence.}
  \item{nblocks}{Number of blocks.}
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
assessStationarity(sequence,1)
divergenceTest(sequence, mcFit$estimate@transitionMatrix)
}

