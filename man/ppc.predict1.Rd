\name{ppc.predict1}
\alias{ppc.predict1}
\title{ Function to do a   single test set prediction for the  PPC method}
\description{
This function  does test set prediction for PPC method.
It predicts  outcome classes  for a   list of peaks  from a single  test set spectrum.}

\usage{
ppc.predict1(centroid.fit, split.fit, logmz, peaklist.te,  threshold, metric = c("binomial", "euclidean", "absolute"), summ = c("mean", "median"))
}
\arguments{
  \item{centroid.fit}{ Result of a call to ppc.make.centroid.list}
  \item{split.fit}{ Result of a call to ppc.find.splits}
  \item{logmz}{log of m/z values from training data} 
  \item{peaklist.te}{List of peaks from test set-  each component is a matrix of
log m/z values and peak intensities}
  \item{threshold}{ Threshold values to use}
  \item{metric}{ "binomial","euclidean", or  "absolute"}
  \item{summ}{ "mean" or median"}
}
\value{
  \item{yhat }{Matrix of predicted classes}
  \item{threshold }{Threshold values used.}
  \item{numsites }{Number of  sites surviving the threshold  for each shrinkage value}
  \item{sites }{List of sites surviving the threshold  for each shrinkage value}
  \item{ind }{Indicator matrix of event ( peak intensity at site > cutpoint}
  \item{ind0 }{Indicator matrix of event ( peak intensity at site > cutpoint}
  \item{prob }{Matrix of estimated class probabilities}
  \item{ht }{Matrix of peak intensities}
}
\author{ Balasubramanian Narasimhan and Rob Tibshirani }


\examples{
## for a complete worked example of this function in a PPC analysis  see
##  http://www-stat.stanford.edu/~tibs/PPC/Rdist/Rscript.rawdata
}

\keyword{ }% at least one, from doc/KEYWORDS

