\name{ppc.predict.peaks1}
\alias{ppc.predict.peaks1}
\title{ A function to find centroid peaks, in a list of individual peaks from one spectrum }
\description{
 Takes centroid.fit (result of a call to make.centroids.list)
  and looks for these m peaks  in a  peaklist from one spectrum.
   returns   ind - an m-vector  of TRUE/FALSE values
     and ht, the vector of  corresponding peak heights. A peak
is considered present if it is within user.parms\$peak.gap units from the centroid.
}
\usage{
ppc.predict.peaks1(centroid.fit, logmz, peaklist.new)
}

\arguments{
  \item{centroid.fit}{ Result of call to make.centroid.list }
  \item{logmz}{Log of m/z values}
\item{ peaklist.new}{Matrix of peaks from a single spectrum. Rows are
(log m/z value, peak intensity)}
}

\value{
  \item{ind }{indicator vector of presence/absence of peak}
  \item{ht }{Vector of  heights. Note: peak is not present if ind=0, even if ht is >0}
}
\author{ Balasubramanian Narasimhan and Rob Tibshirani }


\examples{
## for a complete worked example of this function in a PPC analysis  see
##  http://www-stat.stanford.edu/~tibs/PPC/Rdist/Rscript.rawdata
}

\keyword{ }% at least one, from doc/KEYWORDS


