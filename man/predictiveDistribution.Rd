\name{predictiveDistribution}
\alias{predictiveDistribution}
\title{
Function to compute the probability of observing a new data set, given a data set
}
\description{
The function computes the probability of observing a new data set using information from a given data set. Additionally, hyperparameters can be provided.
}
\usage{
predictiveDistribution(stringchar, newData, hyperparam = matrix(1, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stringchar}{
This is the data using which the Bayesian inference is performed.
}
  \item{newData}{
This is the data whose predictive probability is computed.  
}
  \item{hyperparam}{
This determines the shape of the prior distribution of the parameters.
}
}
\details{
The underlying method is Bayesian inference. The probability is computed by averaging the likelihood of the new data with respect to the posterior. Since the method assumes conjugate priors, the result can be represented in a closed form (see the vignette for more details), which is what is returned.
}
\value{
A single floating point number is returned which is the probability of observing the new data.
}
\references{
Inferring Markov Chains: Bayesian Estimation, Model Comparison, Entropy Rate, and Out-of-Class Modeling, Christopher C. Strelioff, James P. Crutchfield, Alfred Hubler, Santa Fe Institute
}
\author{
Sai Bhargav Yalamanchi
}

\seealso{
\code{\link{markovchainFit}}
}
\examples{
sequence<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", 
"b", "b", "b", "a")
predProb <- predictiveDistribution(sequence[1:10], sequence[11:17], hyperparam = matrix(1, 1))
}