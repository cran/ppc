\name{ppc.peaks}
\alias{ppc.peaks}
\title{Find local maxima}
\description{ Finds the local maxima in a vector.
    This is an R  implementation of the Splus function peaks.
. Note that the span parameter is a proportion between 0 and 1, rather than 
the number of x values (as in the Splus function).
Note also that it only handles a vector as input.
}
\usage{
ppc.peaks(x, span)
}
\arguments{
  \item{x}{ A vector. Peaks will find the local
          maxima in x. }
  \item{span}{A peak is defined as an element in a sequence which is greater
          than all other elements within a window of width length(x)*span centered
          at that element.
  }
}
\details{
 All elements within a halfspan of the end of a sequence or within a
   halfspan of a missing value are FALSE.
}
\value{
vector of logical values, indicating whether there is a peak at each location
}
\author{ Balasubramanian Narasimhan and Rob Tibshirani }

\examples{
x<-rnorm(1000)
a<-ppc.peaks(x, .02)
}

\keyword{ }% at least one, from doc/KEYWORDS


