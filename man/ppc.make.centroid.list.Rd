\name{ppc.make.centroid.list}
\alias{ppc.make.centroid.list}
\title{Function to make a list of peak centroids, from a set of peaks from different spectra}
\description{
This function starts with a list of peaks from a collection of spectra, and does a one
dimensional hierarchical clustering. It  then cuts off the dendogram at height
user.parms\$peak.gap, forming clusters of peaks. The medoids of each cluster form
hte final list of peak centroids.
}
\usage{
ppc.make.centroid.list(data, user.parms)
}
\arguments{
  \item{data}{ List containing mass spc data }
  \item{user.parms}{List of user-defined parameters}
}
\value{
  \item{cent }{Matrix of centroids, one per row}
  \item{peaklist }{Peaklist used}
  \item{clust.tree }{Dendrogram from hclust.1d}
  \item{all.peaks }{Vector of m/z values of all peaks used}
  \item{peak.gap }{Peak width used}
  \item{recluster}{Not currently used}
}
\author{ Balasubramanian Narasimhan and Rob Tibshirani }

\examples{
## for a complete worked example of this function in a PPC analysis  see
##  http://www-stat.stanford.edu/~tibs/PPC/Rdist/Rscript.rawdata
}

\keyword{ }% at least one, from doc/KEYWORDS

