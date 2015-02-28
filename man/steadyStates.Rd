\name{steadyStates}
\alias{steadyStates}
\title{
Stationary states of a \code{markovchain} objeect
}
\description{
This method returns the stationary vector in matricial form of a markovchain object.
}
\usage{
steadyStates(object)
}
\arguments{
  \item{object}{A discrete \code{markovchain} object}
}
\value{A matrix corresponding to the stationary states}
\references{
A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall 2010
}
\author{
Giorgio Spedicato
}
\seealso{
\code{\linkS4class{markovchain}}
}
\note{
The steady states are identified starting from which eigenvectors correspond to identity eigenvalues and then normalizing them to sum up to unity.
}
\examples{
statesNames=c("a","b","c")
markovB<-new("markovchain", states=statesNames, transitionMatrix=
          matrix(c(0.2,0.5,0.3,
                   0,1,0,
                   0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(statesNames,statesNames)
                 ))
steadyStates(markovB)
}

