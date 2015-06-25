\name{priorDistribution}
\alias{priorDistribution}

\title{
Prior Distribution
}
\description{
Function to evaluate the prior probability of a transition matrix. It is based on conjugate priors and therefore a Dirichlet distribution is used to model the transitions of each state.
}
\usage{
priorDistribution(transMatr, hyperparam = matrix())
}
\arguments{
  \item{transMatr}{
The transition matrix whose probability is the parameter of interest.
}
  \item{hyperparam}{
The hyperparam matrix (optional). If not provided, a default value of 1 is assumed for each and therefore the resulting probability distribution is uniform.
}
}
\details{
The states (dimnames) of the transition matrix and the hyperparam may be in any order. 
}
\value{
The log of the probabilities for each state is returned in a numeric vector. Each number in the vector represents the probability (log) of having a probability transition vector as specified in corresponding the row of the transition matrix.
}
\references{
Yalamanchi SB, Spedicato GA (2015). Bayesian Inference of First Order Markov Chains. R
package version 0.2.5
}
\author{
Sai Bhargav Yalamanchi, Giorgio Spedicato
}
\note{
This function can be used in conjunction with inferHyperparam. For example, if the user has a prior data set and a prior transition matrix, he can infer the hyperparameters using inferHyperparam and then compute the probability of their prior matrix using the inferred hyperparameters with priorDistribution.
}


\seealso{
\code{\link{predictiveDistribution}}, \code{\link{inferHyperparam}}
}
\examples{
priorDistribution(matrix(c(0.5, 0.5, 0.5, 0.5), 
                        nrow = 2, 
                        dimnames = list(c("a", "b"), c("a", "b"))), 
                 matrix(c(2, 2, 2, 2), 
                        nrow = 2, 
                        dimnames = list(c("a", "b"), c("a", "b"))))
}