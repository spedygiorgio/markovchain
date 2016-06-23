\name{markovchainList-class}
\Rdversion{1.1}
\docType{class}
\alias{markovchainList-class}
\alias{[[,markovchainList-method}
\alias{dim,markovchainList-method}
\alias{predict,markovchainList-method}
\alias{print,markovchainList-method}
\alias{show,markovchainList-method}
\title{Class \code{"markovchainList"}}
\description{
A class to handle non - homogeneous Markov chains
}
\section{Objects from the Class}{
A \code{markovchainlist} is a list of \code{markovchain} objects. They can be used to model non - homogeneous discrete time Markov Chains, when transition probabilities (and possible states) change by time.
}
\section{Slots}{
  \describe{
    \item{\code{markovchains}:}{Object of class \code{"list"}: a 
    list of markovchains}
    \item{\code{name}:}{Object of class \code{"character"}: optional 
    name of the class }
  }
}
\section{Methods}{
  \describe{
    \item{[[}{\code{signature(x = "markovchainList")}: extract the
    i-th \code{markovchain} }
    \item{dim}{\code{signature(x = "markovchainList")}: number 
    of \code{markovchain} underlying the matrix }
	  \item{predict}{\code{signature(object = "markovchainList")}: predict 
	  from a \code{markovchainList} }
    \item{print}{\code{signature(x = "markovchainList")}: prints the list of markovchains }
    \item{show}{\code{signature(object = "markovchainList")}: same as \code{print} }
   }
}
\references{
A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall 2010
}
\author{
Giorgio Spedicato
}
\note{
The class consists in a list of \code{markovchain} objects. 
It can help to deal with non - homogeneous Markov chains.
}
\seealso{
\code{\linkS4class{markovchain}}
}
\examples{
showClass("markovchainList")
#define a markovchainList
statesNames=c("a","b")

mcA<-new("markovchain",name="MCA", transitionMatrix=matrix(c(0.7,0.3,0.1,0.9),
byrow=TRUE, nrow=2, dimnames=list(statesNames,statesNames)))

mcB<-new("markovchain", states=c("a","b","c"), name="MCB",
transitionMatrix=matrix(c(0.2,0.5,0.3,0,1,0,0.1,0.8,0.1),
nrow=3, byrow=TRUE))

mcC<-new("markovchain", states=c("a","b","c","d"), name="MCC",
  	transitionMatrix=matrix(c(0.25,0.75,0,0,0.4,0.6,
    0,0,0,0,0.1,0.9,0,0,0.7,0.3), nrow=4, byrow=TRUE)
)
mcList<-new("markovchainList",markovchains=list(mcA, mcB, mcC), 
name="Non - homogeneous Markov Chain")
}
\keyword{classes}
