\name{preproglucacon}
\alias{preproglucacon}
\docType{data}
\title{
Preprogluccacon DNA protein bases sequences
}
\description{
Sequence of bases for preproglucacon DNA protein
}
\usage{data(preproglucacon)}
\format{
  A data frame with 1572 observations on the following 2 variables.
  \describe{
    \item{\code{V1}}{a numeric vector, showing original coding}
    \item{\code{preproglucacon}}{a character vector, showing initial of DNA bases (Adenine, Cytosine, Guanine, Thymine)}
  }
}
\source{
Avery Henderson
}
\references{
Averuy Henderson, Fitting markov chain models on discrete time series such as DNA sequences
}
\examples{
data(preproglucacon)
preproglucaconMc<-markovchainFit(data=preproglucacon$preproglucacon)
}
\keyword{datasets}


