\name{ppc.cv}
\alias{ppc.cv}
\title{ Cross-validation for PPC analysis }
\description{
This function does K-fold cross-validation for PPC analysis 
}
\usage{
ppc.cv(ppc.fit, data, user.parms)
}
\arguments{
  \item{ppc.fit}{ Result of call to ppc.predict }
  \item{data}{List containing mass spec data}
  \item{user.parms}{ List of user defiend parameters}
}
\details{
}
\value{
\item{err}{CV error rate for each threshold value }
\item{se}{se of CV error rate}
\item{confusion}{ Cnfusion matrix for each threshold value}
\item{threshold}{Threshold vlaues used}
\item{yhat}{Predicted values from CV}
\item{prob}{Cv probabilities}
\item{y}{Training set outcome values}
\item{folds}{Indices defining CV folds}
\item{numsites}{Number of m/z sites surviving shrinkage  at each threshold value}
}
\author{ Balasubramanian Narasimhan and Rob Tibshirani }


\examples{
## for a complete worked example of this function in a PPC analysis  see
##  http://www-stat.stanford.edu/~tibs/PPC/Rdist/Rscript.rawdata
}

\keyword{ }% at least one, from doc/KEYWORDS


