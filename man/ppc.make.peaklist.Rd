\name{ppc.make.peaklist}
\alias{ppc.make.peaklist}
\title{Function to extract peaks from raw mass spec data }
\description{
This  function to extracts peaks from raw mass spec data. It uses
a very simple peak finder "peaks"- looking  for sites where the intensity is higher than
it is elsewhere in a user-defined window surrounding that site.
}
\usage{
ppc.make.peaklist(data, user.parms)
}
\arguments{
  \item{data}{List containing raw mass spec data}
  \item{user.parms}{ List of user parameters}
}
\value{
List of peaks- one component per spectrum. Each component is a matrix of log m/z values
and peak intensity.
}
\author{ Balasubramanian Narasimhan and Rob Tibshirani }

\examples{
## for a complete worked example of this function in a PPC analysis  see
##  http://www-stat.stanford.edu/~tibs/PPC/Rdist/Rscript.rawdata
}

\keyword{ }% at least one, from doc/KEYWORDS

