\name{inferHyperparam}
\alias{inferHyperparam}

\title{
Function to infer the hyperparameters for Bayesian inference from an a priori matrix or a data set
}
\description{
Since the Bayesian inference approach implemented in the package is based on conjugate priors, hyperparameters must be provided to model the prior probability distribution of the chain parameters. The hyperparameters are inferred from a given a priori matrix under the assumption that the matrix provided corresponds to the mean (expected) values of the chain parameters. A scaling factor vector must be provided too. Alternatively, the hyperparameters can be inferred from a data set. 
}
\usage{
inferHyperparam(transMatr = matrix(), scale = numeric(), data = character())
}

\arguments{
  \item{transMatr}{
A valid transition matrix, with dimension names.
}
  \item{scale}{
A vector of scaling factors, each element corresponds to the row names of the provided transition matrix transMatr, in the same order. 
}
  \item{data}{
A data set from which the hyperparameters are inferred.  
}
}
\details{
transMatr and scale need not be provided if data is provided.
}
\value{
Returns the hyperparameter matrix in a list.
}
\note{
The hyperparameter matrix returned is such that the row and column names are sorted alphanumerically, and the elements in the matrix are correspondingly permuted. 
}

\author{
Sai Bhargav Yalamanchi, Giorgio Spedicato
}

\seealso{
\code{\link{markovchainFit}}, \code{\link{predictiveDistribution}}
}

\examples{
data(rain, package = "markovchain")
inferHyperparam(data = rain$rain)

weatherStates <- c("sunny", "cloudy", "rain")
weatherMatrix <- matrix(data = c(0.7, 0.2, 0.1, 
                                 0.3, 0.4, 0.3, 
                                 0.2, 0.4, 0.4), 
                        byrow = TRUE, nrow = 3, 
                        dimnames = list(weatherStates, weatherStates))
inferHyperparam(transMatr = weatherMatrix, scale = c(10, 10, 10))
}