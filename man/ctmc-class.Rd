\name{ctmc-class}
\Rdversion{1.1}
\docType{class}
\alias{ctmc-class}
\alias{dim,ctmc-method}
\alias{initialize,ctmc-method}
\alias{states,ctmc-method}
\alias{steadyStates,ctmc-method}
\alias{plot,ctmc,missing-method}

\title{Class \code{"ctmc"}}
\description{
The S4 class that describes \code{ctmc} (continuous time Markov chain) objects.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ctmc", states, byrow, generator, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{states}:}{Name of the states. Must be the same of \code{colnames} and \code{rownames} of the generator matrix}
    \item{\code{byrow}:}{Binary flag.}
    \item{\code{generator}:}{Square generator matrix}
    \item{\code{name}:}{Optional character name of the Markov chain}
  }
}
\section{Methods}{
  \describe{
    \item{dim}{\code{signature(x = "ctmc")}: method to get the size}
    \item{initialize}{\code{signature(.Object = "ctmc")}: initialize method }
    \item{states}{\code{signature(object = "ctmc")}: states method. }
    \item{steadyStates}{\code{signature(object = "ctmc")}: method to get the steady state vector. }
    \item{plot}{\code{signature(x = "ctmc", y = "missing")}: plot method for \code{ctmc} objects }
	 }
}
\references{
Introduction to Stochastic Processes with Applications in the Biosciences (2013), David F. Anderson, University of Wisconsin at Madison
}
\author{Sai Bhargav Yalamanchi, Giorgio Spedicato}
\note{
\code{ctmc} objects are written using S4 Classes.
}
\section{Warning }{Validation method is used to assess whether either columns or rows totals to zero. Rounding is used up to 5th decimal. If state names are not properly defined for a generator  \code{matrix}, coercing to \code{ctmc} object leads to overriding states name with artificial "s1", "s2", ... sequence }
\seealso{
\code{\link{generatorToTransitionMatrix}},\code{\link{rctmc}}
}
\examples{
energyStates <- c("sigma", "sigma_star")
byRow <- TRUE
gen <- matrix(data = c(-3, 3,
                       1, -1), nrow = 2,
              byrow = byRow, dimnames = list(energyStates, energyStates))
molecularCTMC <- new("ctmc", states = energyStates, 
                 byrow = byRow, generator = gen, 
                 name = "Molecular Transition Model")

steadyStates(molecularCTMC)
\dontrun{plot(molecularCTMC)}
}
\keyword{classes}
