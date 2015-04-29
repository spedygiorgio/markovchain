\name{blanden}
\alias{blanden}
\docType{data}
\title{
Mobility between income quartiles
}
\description{
This table show mobility between income quartiles for father and sons for the 1970 cohort born
}
\usage{data(blanden)}
\format{
  The format is:
 table [1:4, 1:4] 0.38 0.25 0.21 0.16 0.29 0.28 0.26 0.17 0.22 0.26 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:4] "Bottom" "2nd" "3rd" "Top"
  ..$ : chr [1:4] "Bottom" "2nd" "3rd" "Top"
}
\details{
The rows represent father's income quartile when the son is aged 16, whilst the columns represent sons' income quartiles when he is aged 30 (in 2000).
}
\source{
Personal elaborations from references
}
\references{
Jo Blanden, Paul Gregg and Stephen Machin, Intergenerational Mobility in Europe and North America, Center for Economic Performances (2005)
}
\examples{
data(blanden)
mobilityMc<-as(blanden, "markovchain")
}
\keyword{datasets}
