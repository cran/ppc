\name{ppc.find.splits}
\alias{ppc.find.splits}
\title{ Function to find best discriminating split points for training data in mass spec}
\description{
 Find best discriminating split points for training data (to separate the classes).
  Takes centroids.fit- result of call to make.centroids.list
   and peaks.fit- result of call to predict.peaks
}
\usage{
ppc.find.splits(centroid.fit, peak.fit, data, user.parms)
}
\arguments{
  \item{centroid.fit}{ Result of call to ppc.make.centroid.list }
  \item{peak.fit}{  Result of call to  ppc.predict.peaks }
  \item{data}{ List containing mass spec data}
  \item{user.parms}{ List of user defined parameters}
}

\value{
\item{prhat}{Proportion of samples beyond optimal cutpoint  in each outcome class}
\item{pr}{Proportion of samples beyond cutpoints in each outcome class}
\item{n.class}{number of samples in each outcome class}
\item{cutpoints}{Cutpoints (split points) tried}
\item{cuthat}{Optimal cut points}
\item{prclose}{Indicators for split points with prob difference within
  10 percent of that of the optimal split point}
\item{nsplits}{Number of cutpoints tried}
\item{fix.at.one}{Was the optimal cutpoint fixed at one? (i.e no peak vs peak)}
}
\author{ Balasubramanian Narasimhan and Rob Tibshirani }

\examples{
## for a complete worked example of this function in a PPC analysis  see
##  http://www-stat.stanford.edu/~tibs/PPC/Rdist/Rscript.rawdata
}

\keyword{ }% at least one, from doc/KEYWORDS


