\name{markovchain-class}
\Rdversion{1.1}
\docType{class}
\alias{markovchain-class}
\alias{*,markovchain,markovchain-method}
\alias{*,markovchain,matrix-method}
\alias{*,markovchain,numeric-method}
\alias{*,matrix,markovchain-method}
\alias{*,numeric,markovchain-method}
\alias{[,markovchain,ANY,ANY,ANY-method}
\alias{^,markovchain,numeric-method}
\alias{==,markovchain,markovchain-method}
\alias{!=,markovchain,markovchain-method}
\alias{absorbingStates,markovchain-method}
\alias{conditionalDistribution,markovchain-method}
\alias{canonicForm,markovchain-method}
\alias{coerce,data.frame,markovchain-method}
\alias{coerce,markovchain,data.frame-method}
\alias{coerce,table,markovchain-method}
\alias{coerce,markovchain,igraph-method}
\alias{coerce,markovchain,matrix-method}
\alias{coerce,markovchain,sparseMatrix-method}
\alias{coerce,sparseMatrix,markovchain-method}
\alias{coerce,matrix,markovchain-method}
\alias{coerce,msm,markovchain-method}
\alias{coerce,msm.est,markovchain-method}
\alias{coerce,etm,markovchain-method}
\alias{dim,markovchain-method}
\alias{initialize,markovchain-method}
\alias{names,markovchain-method}
\alias{names<-,markovchain-method}
\alias{plot,markovchain,missing-method}
\alias{predict,markovchain-method}
\alias{print,markovchain-method}
\alias{show,markovchain-method}
\alias{summary,markovchain-method}
\alias{sort,markovchain-method}
\alias{t,markovchain-method}
\title{Class \code{"markovchain"}}
\description{
The S4 class that describes \code{markovchain} objects.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("markovchain", states, byrow, transitionMatrix, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{states}:}{Name of the states. Must be the same of \code{colnames} and \code{rownames} of the transition matrix}
    \item{\code{byrow}:}{Binary flag.}
    \item{\code{transitionMatrix}:}{Square transition matrix}
    \item{\code{name}:}{Optional character name of the Markov chain}
  }
}
\section{Methods}{
  \describe{
    \item{*}{\code{signature(e1 = "markovchain", e2 = "markovchain")}: multiply two \code{markovchain} objects}
    \item{*}{\code{signature(e1 = "markovchain", e2 = "matrix")}: markovchain by matrix multiplication}
    \item{*}{\code{signature(e1 = "markovchain", e2 = "numeric")}: markovchain by numeric vector multiplication }
    \item{*}{\code{signature(e1 = "matrix", e2 = "markovchain")}: matrix by markov chain}
    \item{*}{\code{signature(e1 = "numeric", e2 = "markovchain")}: numeric vector by \code{markovchain} multiplication   }
    \item{[}{\code{signature(x = "markovchain", i = "ANY", j = "ANY", drop = "ANY")}: ... }
    \item{^}{\code{signature(e1 = "markovchain", e2 = "numeric")}: power of a \code{markovchain} object}
    \item{==}{\code{signature(e1 = "markovchain", e2 = "markovchain")}: equality of two \code{markovchain} object}
    \item{!=}{\code{signature(e1 = "markovchain", e2 = "markovchain")}: non-equality of two \code{markovchain} object}

    \item{absorbingStates}{\code{signature(object = "markovchain")}: method to get absorbing states }
	\item{canonicForm}{\code{signature(object = "markovchain")}: return a \code{markovchain} object into canonic form }
    \item{coerce}{\code{signature(from = "markovchain", to = "data.frame")}: coerce method from markovchain to \code{data.frame}}
    \item{conditionalDistribution}{\code{signature(object = "markovchain")}: returns the conditional probability of subsequent states given a state}
    \item{coerce}{\code{signature(from = "data.frame", to = "markovchain")}: coerce method from \code{data.frame} to \code{markovchain}}
    \item{coerce}{\code{signature(from = "table", to = "markovchain")}: coerce method from \code{table} to \code{markovchain} }
    \item{coerce}{\code{signature(from = "msm", to = "markovchain")}: coerce method from \code{msm} to \code{markovchain} }
    \item{coerce}{\code{signature(from = "msm.est", to = "markovchain")}: coerce method from \code{msm.est} (but only from a Probability Matrix) to \code{markovchain} }
    \item{coerce}{\code{signature(from = "etm", to = "markovchain")}: coerce method from \code{etm} to \code{markovchain} }
    \item{coerce}{\code{signature(from = "sparseMatrix", to = "markovchain")}: coerce method from \code{sparseMatrix} to \code{markovchain} }
    \item{coerce}{\code{signature(from = "markovchain", to = "igraph")}: coercing to \code{igraph} objects }
    \item{coerce}{\code{signature(from = "markovchain", to = "matrix")}: coercing to \code{matrix} objects }
    \item{coerce}{\code{signature(from = "markovchain", to = "sparseMatrix")}: coercing to \code{sparseMatrix} objects }
    \item{coerce}{\code{signature(from = "matrix", to = "markovchain")}: coercing to \code{markovchain} objects from \code{matrix} one }
    \item{dim}{\code{signature(x = "markovchain")}: method to get the size}
    \item{names}{\code{signature(x = "markovchain")}: method to get the names of states}
    \item{names<-}{\code{signature(x = "markovchain", value = "character")}: method to set the names of states}
    \item{initialize}{\code{signature(.Object = "markovchain")}: initialize method }
    \item{plot}{\code{signature(x = "markovchain", y = "missing")}: plot method for \code{markovchain} objects }
	\item{predict}{\code{signature(object = "markovchain")}: predict method }
    \item{print}{\code{signature(x = "markovchain")}: print method. }
    \item{show}{\code{signature(object = "markovchain")}: show method. }
    \item{sort}{\code{signature(x = "markovchain", decreasing=FALSE)}: sorting the transition matrix. }
    \item{states}{\code{signature(object = "markovchain")}: returns the names of states (as \code{names}. }
    \item{steadyStates}{\code{signature(object = "markovchain")}: method to get the steady vector. }
	 \item{summary}{\code{signature(object = "markovchain")}: method to summarize structure of the markov chain }
	\item{transientStates}{\code{signature(object = "markovchain")}: method to get the transient states. }
    \item{t}{\code{signature(x = "markovchain")}: transpose matrix }
    \item{transitionProbability}{\code{signature(object = "markovchain")}: transition probability }
	 }
}
\references{
A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall 2010
}
\author{Giorgio Spedicato}
\note{
\code{markovchain} object are written in S4 Classes.
}
\section{Warning }{Validation method is used to assess whether either columns or rows totals to one. Rounding is used up to 5th decimal. If state names are not properly defined for a probability  \code{matrix}, coercing to \code{markovhcain} object leads to overriding states name with artificial "s1", "s2", ... sequence. In addition, operator overloading has been applied for \eqn{+,*,^,==,!=} operators. }
\seealso{
\code{\link{markovchainSequence}},\code{\link{markovchainFit}}
}
\examples{
#show markovchain definition
showClass("markovchain")
#create a simple Markov chain
transMatr<-matrix(c(0.4,0.6,.3,.7),nrow=2,byrow=TRUE)
simpleMc<-new("markovchain", states=c("a","b"),
transitionMatrix=transMatr, 
name="simpleMc")
#power
simpleMc^4
#some methods
steadyStates(simpleMc)
absorbingStates(simpleMc)
simpleMc[2,1]
t(simpleMc)
is.irreducible(simpleMc)
#conditional distributions
conditionalDistribution(simpleMc, "b")
#example for predict method
sequence<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", "b", "b", "b", "a")
mcFit<-markovchainFit(data=sequence)
predict(mcFit$estimate, newdata="b",n.ahead=3)
#direct conversion
myMc<-as(transMatr, "markovchain")

#example of summary
summary(simpleMc)
\dontrun{plot(simpleMc)}
}
\keyword{classes}
