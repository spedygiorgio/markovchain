\name{generatorToTransitionMatrix}
\alias{generatorToTransitionMatrix}

\title{
Function to obtain the transition matrix from the generator.
}
\description{
The transition matrix of the embedded DTMC is inferred from the CTMC's generator.
}
\usage{
generatorToTransitionMatrix(gen, byrow = TRUE)
}

\arguments{
  \item{gen}{
The generator matrix.
}
  \item{byrow}{
Flag to determine if rows (columns) sum to 0.
}
}

\value{
Returns the transition matrix.
}
\references{
Introduction to Stochastic Processes with Applications in the Biosciences (2013), David F. Anderson, University of Wisconsin at Madison
}
\author{
Sai Bhargav Yalamanchi
}

\seealso{
\code{\link{rctmc}},\code{\link{ctmc-class}}
}
\examples{
energyStates <- c("sigma", "sigma_star")
byRow <- TRUE
gen <- matrix(data = c(-3, 3,
                       1, -1), nrow = 2,
              byrow = byRow, dimnames = list(energyStates, energyStates))
generatorToTransitionMatrix(gen)
}

