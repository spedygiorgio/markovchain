\name{multinomialConfidenceItervals}
\alias{multinomialConfidenceItervals}
\title{
A function to compute multinomial confidence intervals of DTMC
}
\description{
These functions return multinomial confidence intervals of the \code{markovchain} objects.
}
\usage{
multinomialConfidenceItervals(transitionMatrix, countsTransitionMatrix, confidencelevel)
}
\arguments{
  \item{transitionMatrix}{Transition matrix.}
  \item{countsTransitionMatrix}{Counts transition matrix.}
  \item{confidencelevel}{Confidence level.}
}
\value{A list with two matrices}
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
seq<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
mcfit<-markovchainFit(data=seq,byrow=TRUE)
seqmat<-createSequenceMatrix(seq)
multinomialConfidenceItervals(mcfit$estimate@transitionMatrix, seqmat, 0.95)
}

