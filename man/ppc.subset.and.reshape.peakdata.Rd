\name{ppc.subset.and.reshape.peakdata}
\alias{ppc.subset.and.reshape.peakdata}
\title{Function to subset and reshape  mass spec peak data }
\description{
This function subsets   mass spec peak data, using the m/z limits mz.min and mz.max in
user.parms. If there are batches, it  also reshapes the data, concatenating the
spectra in different batches for a given patient. Thus for each patient it produces one 
long vector of values. The m/z values are strung out as well. Using this trick
we can essentially ignore the  presence of batches for the rest of the analysis,
}
\usage{
ppc.subset.and.reshape.peakdata(data, user.parms)
}
\arguments{
  \item{data}{ List containing the mass spec data }
  \item{user.parms}{ List of user parameters}
}
\value{
List containing the subset and reshaped mass spec data
}
\author{ Balasubramanian Narasimhan and Rob Tibshirani }


\examples{
## for a complete worked example of this function in a PPC analysis  see
##  http://www-stat.stanford.edu/~tibs/PPC/Rdist/Rscript.peakdata
}

\keyword{ }% at least one, from doc/KEYWORDS

