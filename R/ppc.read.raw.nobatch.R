##
## Read raw data from a directory and return a list containing
## -- a matrix of x
## -- a vector of the mz values
## -- the list of file names
##

ppc.read.raw.nobatch <- function(directory, mz = NULL) {
  datafiles.list <- ppc.xl.get.names.of.files(directory)
  
  ##
  ## First determine the dimensions of the data matrix
  ##

  if (is.null(mz)) {
    pat1 <- read.table(datafiles.list[1])
    mz <- pat1[,1]
  }
  
  xtr <- matrix(NA, nrow=length(mz), ncol=length(datafiles.list))

  ##
  ## Now read in all the data
  ##
  for (j in 1:length(datafiles.list)) {
    if (ppc.options$debug) print(paste("Reading file", datafiles.list[j]))
    
    temp <-  read.table(datafiles.list[j], sep = ",")
    
    xtr[, j] <- approx(temp[, 1], temp[, 2], xout = mz)$y
  }
  
  return(list(xtr = xtr,
              mz = mz,
              filenames=datafiles.list))
}
