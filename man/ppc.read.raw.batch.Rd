\name{ppc.read.raw.batch}
\alias{ppc.read.raw.batch}
\title{A function to read in raw protein mass spec data, with  batches }
\description{This function  reads in raw protein mass spec data  from a directory.
The directory is assumed to have on subdirectory per class (eg control or disease).
The subdirectories like "control" have further subdirectories, one per batch.
So the structure looks lik control/batch1/file.csv, control/batch2/file.csv,
disease/batch1/file.csv, etc.
There is one  comma-separated (csv) file.csv in the subdirectory per spectrum,  
having lines of the form m/ value, intensity (one line per m/z site)}
\usage{
ppc.read.raw.batch(dir, batches, mz = NULL)
}
\arguments{
  \item{dir}{ Name of directory containing the data }
  \item{batches}{Vector of character names of batches }
  \item{mz}{Optional vector of m/z values. Default NULL. If NULL, m/z values
are read in from files. Otherwise the  values  in mz are used.
   }
}
\value{
\item{xtr}{Matrix of intensities- one row per m/z site, one col per spectrum}
  \item{mz}{m/z values }
\item{filenames}{List of filenames read in}
}
\author{ Balasubramanian Narasimhan and Rob Tibshirani }


\examples{
## for a complete worked example of this function in a PPC analysis  see
##  http://www-stat.stanford.edu/~tibs/PPC/Rdist/Rscript.rawdata
}

\keyword{ }% at least one, from doc/KEYWORDS
