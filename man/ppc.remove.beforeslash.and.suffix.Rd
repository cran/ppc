\name{ppc.remove.beforeslash.and.suffix}
\alias{ppc.remove.beforeslash.and.suffix}
\title{ Remove characters efore a slash and at end  }
\description{
This function takes a character file name or vector of character file names,
and removes the character before the rightmost slash, and the suffix after the
last dot. Eg "file/foo/junk.csv" becomes "junk"
}
\usage{
ppc.remove.beforeslash.and.suffix(x)
}
\arguments{
  \item{x}{ Character string or vector of character strings. }
}
\details{
}
\value{
Character string or vector of character strings, modified in the manner descibed
above.
}
\author{ Balasubramanian Narasimhan and Rob Tibshirani }


\examples{
## for a complete worked example of this function in a PPC analysis  see
##  http://www-stat.stanford.edu/~tibs/PPC/Rdist/Rscript.rawdata
}

\keyword{ }% at least one, from doc/KEYWORDS
