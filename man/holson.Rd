\name{holson}
\alias{holson}
\docType{data}
\title{
Holson data set
}
\description{
A data set containing 1000 life histories trajectories and a categorical status (1,2,3) observed on eleven 
evenly spaced steps.
}
\usage{data("holson")}
\format{
  A data frame with 1000 observations on the following 12 variables.
  \describe{
    \item{\code{id}}{unique id}
    \item{\code{time1}}{observed status at i-th time}
    \item{\code{time2}}{observed status at i-th time}
    \item{\code{time3}}{observed status at i-th time}
    \item{\code{time4}}{observed status at i-th time}
    \item{\code{time5}}{observed status at i-th time}
    \item{\code{time6}}{observed status at i-th time}
    \item{\code{time7}}{observed status at i-th time}
    \item{\code{time8}}{observed status at i-th time}
    \item{\code{time9}}{observed status at i-th time}
    \item{\code{time10}}{observed status at i-th time}
    \item{\code{time11}}{observed status at i-th time}
  }
}
\details{
The example can be used to fit a \code{markovchain} or a \code{markovchainList} object.
}
\source{
Private communications
}
\references{
Private communications
}
\examples{
data(holson)
head(holson)
}
\keyword{datasets}
