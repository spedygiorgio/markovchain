\name{rain}
\alias{rain}
\docType{data}
\title{
Alofi island daily rainfall
}
\description{
Rainfall measured in Alofi Island
}
\usage{data(rain)}
\format{
  A data frame with 1096 observations on the following 2 variables.
  \describe{
    \item{\code{V1}}{a numeric vector, showing original coding}
    \item{\code{rain}}{a character vector, showing daily rainfall millilitres brackets}
  }
}
\source{
Avery Henderson
}
\references{
Avery Henderson, Fitting markov chain models on discrete time series such as DNA sequences
}
\examples{
data(rain)
rainMc<-markovchainFit(data=rain$rain)
}
\keyword{datasets}
