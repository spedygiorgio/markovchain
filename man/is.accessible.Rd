\name{is.accessible}
\alias{is.accessible}
\title{Verify if a state j is reachable from state i.}
\description{
This function verifies if a state is reachable from another, i.e., 
if exists a path that leads to state j leaving from state i with positive probability
}
\usage{
is.accessible(object, from, to)
}
\arguments{
  \item{object}{
A \code{markovchain} object.
}
  \item{from}{
The name of state "i" (beginning state).
}
  \item{to}{
The name of state "j" (ending state).
}
}
\details{
If wraps and internal function named \code{.commStatesFinder}.
}
\value{
A boolean value.
}
\references{
James Montgomery, University of Madison
}
\author{
Giorgio Spedicato
}

\seealso{
\code{is.irreducible}
}
\examples{
statesNames<-c("a","b","c")
markovB<-new("markovchain", states=statesNames, transitionMatrix=
          matrix(c(0.2,0.5,0.3,
                   0,1,0,
                   0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(statesNames,statesNames)
                 ))
is.accessible(markovB,"a","c")
}
