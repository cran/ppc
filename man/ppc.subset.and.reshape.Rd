\name{ppc.subset.and.reshape}
\alias{ppc.subset.and.reshape}
\title{Function to subset and reshape raw mass spec data }
\description{
This function subsets  raw mass spec data, using the m/z limits mz.min and mz.max in
user.parms. If there are batches, it  also reshapes the data, concatenating the
spectra in different batches for a given patient. Thus for each patient it produces one 
long column of values. The m/z values are strung out as well. Using this trick
we can essentially ignore the  presence of batches for the rest of the analysis,
}
\usage{
ppc.subset.and.reshape(data, user.parms)
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
##  http://www-stat.stanford.edu/~tibs/PPC/Rdist/Rscript.rawdata
}

\keyword{ }% at least one, from doc/KEYWORDS

