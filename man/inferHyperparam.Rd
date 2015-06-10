\name{inferHyperparam}
\alias{inferHyperparam}
%- Also NEED an '\alias' for EACH other topic documented here.
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
A vector of scaling factors, each element corresponds to the sum of the hyperparameters in the corresponding row of the hyperparameter matrix. 
}
  \item{data}{
A data set from which the hyperparameters are inferred.  
}
}
\details{
transMatr and scale need not be provided if data is provided.
}
\value{
Returns the hyperparameter matrix in a list. The hyperparameter matrix is such that the row and column names are sorted alphanumerically. 
}

\author{
Sai Bhargav Yalamanchi, Giorgio Spedicato
}

\examples{
data(rain, package = "markovchain")
inferHyperparam(data = rain$rain)
}