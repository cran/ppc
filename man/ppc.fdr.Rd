\name{ppc.fdr}
\alias{ppc.fdr}
\title{ Function to estimate False Discovery rates for peaks in PPC analysis}
\description{
 Estimate False Discovery rates for peaks in FDR analysis, using permutations
of the sample labels}
}
\usage{
ppc.fdr(data, centroid.fit, peak.fit, split.fit, ppc.fit, user.parms)
}
\arguments{
  \item{data}{ List containing mass spec data}
 \item{centroid.fit}{ Result of call to ppc.make.centroid.list }
  \item{peak.fit}{ Result of call to  ppc.predict.peaks}
  \item{split.fit}{ Result of call to ppc.find.splits}
\item{ppc.fit}{Result of call to  ppc.predict}
  \item{user.parms}{ List  of user defined parameters}
}
\details{
}
\value{
  \item{results}{Matrix with columns- threshold used, number of peaks found, FDR}
  \item{pi0}{Esimate of proportion of truly null peaks}
\item{threshold}{Vector of thresholds used}
}
\author{ Balasubramanian Narasimhan and Rob Tibshirani }


\examples{
## for a complete worked example of this function in a PPC analysis  see
##  http://www-stat.stanford.edu/~tibs/PPC/Rdist/Rscript.rawdata
}

\keyword{ }% at least one, from doc/KEYWORDS


