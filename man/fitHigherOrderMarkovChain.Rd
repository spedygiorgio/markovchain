\name{fitHigherOrder}
\alias{fitHigherOrder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to fit a higher order Markov chain
}
\description{
Given a sequence of states arising from a stationary state, 
it fits the underlying Markov chain distribution with higher order. 
}
\usage{
fitHigherOrder(data, order = 2)
}

\arguments{
  \item{data}{
A character list.
}
\item{order}{
Markov chain order
}
}
\value{
A list containing lambda.
}
\references{
Ching, W. K., Huang, X., Ng, M. K., & Siu, T. K. (2013). Higher-order markov chains. In Markov Chains (pp. 141-176). Springer US.

Ching, W. K., Ng, M. K., & Fung, E. S. (2008). Higher-order multivariate Markov chains and their applications. Linear Algebra and its Applications, 428(2), 492-507.

package version 0.2.5
}
\author{
Giorgio Spedicato, Tae Seung Kang
}
\note{
This function is written in Rcpp.
}

%\section{Warning }{"mle" method calls \code{createSequenceMatrix} function using \code{sanitize} parameter set to TRUE.}

\seealso{
\code{\link{markovchain}}
}
\examples{
sequence<-c("a", "a", "b", "b", "a", "c", "b", "a", "b", "c", "a", "b", "c", "a"
            , "b", "c", "a", "b", "a", "b")
fitHigherOrder(sequence)
}


