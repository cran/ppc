\name{ppc.read.peaks.nobatch}
\alias{ppc.read.peaks.nobatch}
\title{ Read in mass spec peak data, without batches}
\description{
A function to read in  protein mass spec peaks  data
}
\usage{
ppc.read.peaks.nobatch(dir)
}
\arguments{
  \item{dir}{ Name of directory containing the data }
}
\details{
}
\value{
\item{peaklist}{List of peaks for each sample. Each component is a matrix- one row per m/z site, consiting of log m/z and peak intensity}
\item{filenames}{List of filenames read in}
  \item{logmz}{log m/z values }
}
\author{ Balasubramanian Narasimhan and Rob Tibshirani }


\examples{
## for a complete worked example of this function in a PPC analysis  see
##  http://www-stat.stanford.edu/~tibs/PPC/Rdist/Rscript.peakdata
}

\keyword{ }% at least one, from doc/KEYWORDS
