\name{ppc.plot.hist}
\alias{ppc.plot.hist}
\title{Plot peak histograms from PPC analysis}
\description{
This function plots the  histograms of peaks from a PPC analysis. They are  laid out in order of discriminatory power, starting at the top left and moving down the leftmost column.  At most 25 sites are plotted.
}
\usage{
ppc.plot.hist(peak.fit, ppc.fit, centroid.fit, split.fit, data, first.site = 1, last.site = 25, title.plot = NULL)
}
\arguments{
\item{peak.fit}{ Result of call to  ppc.predict.peaks}
\item{ppc.fit}{ Result of call to  ppc.predict}
 \item{centroid.fit}{ Result of call to ppc.make.centroid.list }
 \item{split.fit}{ Result of call to ppc.find.splits}
  \item{data}{ List containing mass spec data}
  \item{first.site}{ Integer- first site to plot. Default 1}.
  \item{last.site}{Integer- last site to plot. Default 1.}
  \item{title.plot}{Character title for plot. Default NULL}
}
\details{
}
\value{
}
\author{ Balasubramanian Narasimhan and Rob Tibshirani }


\examples{
## for a complete worked example of this function in a PPC analysis  see
##  http://www-stat.stanford.edu/~tibs/PPC/Rscript.rawdata
}

\keyword{ }% at least one, from doc/KEYWORDS
