ppc.predict.peaks <- function(centroid.fit, data)
{
  ## takes centroid.fit (result of a call to make.centroids.list)
  ##  and looks for these m peaks  in peaklist, a list of length n
  ## returns   ind - an m by n matrix of TRUE/FALSE values
  ##       and ht, the matrix of  corresponding peak heights
  
  
peaklist<-data$peaklist
logmz<- data$logmz

  res2 <- NULL
  ht <- NULL
  for(i in 1:length(peaklist)) {
    if (ppc.options$debug) cat(i)
    
    junk <- ppc.predict.peaks1(centroid.fit, logmz,  peaklist.new = 
                           peaklist[[i]])
    
    res2 <- cbind(res2, junk$ind)
    ht <- cbind(ht, junk$ht)
    
  }
  return(list(ind = res2, ht = ht))
}
