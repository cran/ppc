\name{ppc.peak.summary}
\alias{ppc.peak.summary}
\title{ Produce summary of peaks from PPC analysis }
\description{
This function produces a summary of peaks from a PPC analysis 
}
\usage{
ppc.peak.summary(centroid.fit, peak.fit, data, user.parms, split.fit)
}
\arguments{
  \item{centroid.fit}{ Result of call to ppc.make.centroid.list }
  \item{peak.fit}{ Result of call to  ppc.predict.peaks}
  \item{data}{ List containing mass spec data}
  \item{user.parms}{ List  of user defined parameters}
  \item{split.fit}{ Result of call to ppc.find.splits}
}
\details{
}
\value{
Matrix containing on row per peak. Columns are peak.position,
                               min, max of peak.position,
                               number.of.spectra in which the peak occurs,
                               rank.of.peak, in discriminatory power,
                               split point for peak height,
 proportion of samples in class 1 exceeding split point, proportion of samples in class 2 exceeding split point  etc,  pairwise differences in these probabilities,
                               peak info for each sample}
}
\author{ Balasubramanian Narasimhan and Rob Tibshirani }

\examples{
## for a complete worked example of this function in a PPC analysis  see
##  http://www-stat.stanford.edu/~tibs/PPC/Rdist/Rscript.rawdata
}

\keyword{ }% at least one, from doc/KEYWORDS


