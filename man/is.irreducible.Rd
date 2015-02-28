\name{is.irreducible}
\alias{is.irreducible}
\title{Function to check if a Markov chain is irreducible}
\description{
This function verifies whether a \code{markovchain} object transition matrix is composed by only one communicating class.
}
\usage{
is.irreducible(object)
}
\arguments{
  \item{object}{A \code{markovchain} object}
}
\details{
It is based on \code{.communicatingClasses} internal function.
}
\value{
A boolean values.
}
\references{
Feres, Matlab listings for Markov Chains.
}
\author{
Giorgio Spedicato
}
\seealso{
\code{\link{summary}}
}
\examples{
statesNames=c("a","b")
mcA<-new("markovchain", transitionMatrix=matrix(c(0.7,0.3,0.1,0.9),byrow=TRUE, nrow=2, 
                                                  dimnames=list(statesNames,statesNames)
                                                  ))
is.irreducible(mcA)
}
