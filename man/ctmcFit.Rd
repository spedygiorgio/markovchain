\name{ctmcFit}
\alias{ctmcFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to fit a CTMC
}
\description{
This function fits the underlying CTMC give the state transition data and the transition times using the maximum likelihood method (MLE)
}
\usage{
ctmcFit(data, byrow = TRUE, name = "", confidencelevel = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
It is a list of two elements. The first element is a character vector denoting the states. The second is a numeric vector denoting the corresponding transition times. 
}
  \item{byrow}{
Determines if the output transition probabilities of the underlying embedded DTMC are by row.
}
  \item{name}{
Optional name for the CTMC.
}
  \item{confidencelevel}{
Confidence level for the confidence interval construnction. 
}
}
\details{
Note that in data, there must exist an element wise corresponding between the two elements of the list and that data[[2]][1] is always 0.
}
\value{
It returns a list containing the CTMC object and the confidence intervals.
}
\references{
Continuous Time Markov Chains (vignette), Sai Bhargav Yalamanchi, Giorgio Alfredo Spedicato 2015
}
\author{
Sai Bhargav Yalamanchi
}

\seealso{
\code{\link{rctmc}}
}
\examples{
data <- list(c("a", "b", "c", "a", "b", "a", "c", "b", "c"), c(0, 0.8, 2.1, 2.4, 4, 5, 5.9, 8.2, 9))

ctmcFit(data)
}

