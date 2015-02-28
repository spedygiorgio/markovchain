\name{rmarkovchain}
\alias{rmarkovchain}
\alias{markovchainSequence}
\title{
Function to generate a sequence of states from homogeneous or non-homogeneous Markov chains.
}
\description{
Provided any \code{markovchain} or \code{markovchainList} objects, it returns a sequence of 
states coming from the underlying stationary distribution. 
}
\usage{
rmarkovchain(n, object, ...)
markovchainSequence(n, markovchain, t0 = sample(markovchain@states, 1), 
include.t0 = FALSE)
}

\arguments{
\item{n}{
Sample size
}
\item{object}{
Either a \code{markovchain} or a \code{markovchainList} object.
}
\item{\dots}{
additional parameters passed to the internal sampler
}
\item{markovchain}{
The \code{markovchain} object
}
\item{t0}{
The initial state.
}
\item{include.t0}{
Specify if the initial state shall be used.
}
}
\details{
When an homogeneous process is assumed (\code{markovchain} object) a sequence is 
sampled of size n. When an non - homogeneous process is assumed,
n samples are taken but the process is assumed to last from the begin to the end of the 
non-homogeneous markov process.
}
\value{
Either a character vector or a data frame
}
\references{
A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall 2010
}
\author{
Giorgio Spedicato
}
\note{
Check the type of input
}

\seealso{
\code{\link{markovchainFit}}
}
\examples{
#define the Markov chain
statesNames=c("a","b","c")
mcB<-new("markovchain", states=statesNames, transitionMatrix=matrix(c(0.2,0.5,0.3,
0,0.2,0.8,0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(statesNames,statesNames)
                 ))
#show the sequence
outs<-markovchainSequence(n=100,markovchain=mcB, t0="a")
outs2<-rmarkovchain(n=20, object=mcB)
}
