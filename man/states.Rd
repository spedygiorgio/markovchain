\name{states}
\alias{states}
\title{
Defined states of a transition matrix
}
\description{
This method returns the states of a transition matrix.
}
\usage{
states(object)
}
\arguments{
  \item{object}{A discrete \code{markovchain} object}
}
\value{The character vector corresponding to states slot.}
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
statesNames=c("a","b","c")
markovB<-new("markovchain", states=statesNames, transitionMatrix=
          matrix(c(0.2,0.5,0.3,
                   0,1,0,
                   0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(statesNames,statesNames)
                 ))
states(markovB)
}

