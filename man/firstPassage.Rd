\name{firstPassage}
\alias{firstPassage}
\title{First passage across states}
\description{This function compute the first passage probability in states}
\usage{firstPassage(object, state, n)}
\arguments{
  \item{object}{A \code{markovchain} object}
  \item{state}{Initial state}
  \item{n}{Number of rows on which compute the distribution}
}
\details{Based on Feres' Matlab listings}
\value{
A matrix of size 1:n x number of states showing the probability of the 
first time of passage in states to be exactly the number in the row.
}
\references{Renaldo Feres, Notes for Math 450 Matlab listings for Markov chains}
\author{Giorgio Spedicato}
\seealso{
\code{\link{conditionalDistribution}}
}
\examples{
#create a simple Markov chain
simpleMc<-new("markovchain", states=c("a","b"),
transitionMatrix=matrix(c(0.4,0.6,.3,.7),nrow=2,byrow=TRUE))
firstPassage(simpleMc,"b",20)
}