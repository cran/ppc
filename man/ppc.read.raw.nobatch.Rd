\name{ppc.read.raw.nobatch}
\alias{ppc.read.raw.nobatch}
\title{A function to read in raw protein mass spec data, with no batches }
\description{This function  reads in raw protein mass spec data  from a directory.
The directory is assumed to have on subdirectory per class (eg control or disease).
There is one  comma-separated (csv) file in the subdirectory per spectrum,  
having lines of the form m/ value, intensity (one line per m/z site)}
\usage{
ppc.read.raw.nobatch(directory, mz = NULL)
}
\arguments{
  \item{directory}{ Name of directory containing the data }
  \item{mz}{Optional vector of m/z values. Default NULL. If NULL, m/z values
are read in from files. Otherwise the  values  in mz are used.}
}
\value{
  \item{xtr}{Matrix of intensities- one row per m/z site, one col per spectrum}
  \item{mz}{m/z values }
\item{filenames}{List of names of  files that were read in}
}
\author{ Balasubramanian Narasimhan and Rob Tibshirani }


\examples{
## for a complete worked example of this function in a PPC analysis  see
##  http://www-stat.stanford.edu/~tibs/PPC/Rdist/Rscript.rawdata
}

\keyword{ }% at least one, from doc/KEYWORDS
