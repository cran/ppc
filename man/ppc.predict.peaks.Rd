\name{ppc.predict.peaks}
\alias{ppc.predict.peaks}
\title{ A function to find centroid peaks, in a list of individual peaks }
\description{
 Takes centroid.fit (result of a call to make.centroids.list)
  and looks for these m peaks  in peaklist, a list of length n
   returns   ind - an m by n matrix of TRUE/FALSE values
     and ht, the matrix of  corresponding peak heights. A peak
is considered present if it is within user.parms\$peak.gap units from the centroid.
}
\usage{
ppc.predict.peaks(centroid.fit, data)
}

\arguments{
  \item{centroid.fit}{ Result of call to make.centroid.list }
  \item{data}{ List containing the mass spec peaks data }
}

\value{
  \item{ind }{indicator matrix of presence/absence of peak}
  \item{ht }{matrix of  heights. Note: peak is not present if ind=0, even if ht is >0}
}
\author{ Balasubramanian Narasimhan and Rob Tibshirani }


\examples{
## for a complete worked example of this function in a PPC analysis  see
##  http://www-stat.stanford.edu/~tibs/PPC/Rdist/Rscript.rawdata
}

\keyword{ }% at least one, from doc/KEYWORDS


