\name{craigsendi}
\alias{craigsendi}
\docType{data}
\title{
CD4 cells counts on HIV Infects between zero and six month
}
\description{
This is the table shown in Craig and Sendi paper showing zero and six month CD4 cells count in six brakets
}
\usage{data(craigsendi)}
\format{
  The format is:
 table [1:3, 1:3] 682 154 19 33 64 19 25 47 43
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:3] "0-49" "50-74" "75-UP"
  ..$ : chr [1:3] "0-49" "50-74" "75-UP"
}
\details{
Rows represent counts at the beginning, cols represent counts after six months.
}
\source{
Estimation of the transition matrix of a discrete time Markov chain, Bruce A. Craig and Peter P. Sendi, Health Economics 11, 2002.
}
\references{
See source
}
\examples{
data(craigsendi)
csMc<-as(craigsendi, "markovchain")
steadyStates(csMc)
}
\keyword{datasets}
