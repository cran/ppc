ppc.make.peaklist <- function(data, user.parms){
  xtr <- data$xtr
  mz <- data$mz
  logmz<- data$logmz
  
  pk1 <- vector("list",ncol(xtr))
  ii <- 0
  
  for(i in 1:ncol(xtr)) {
    a <- ppc.find.peaks(mz,xtr[,i],user.parms)
    aa<-match(a[,1],mz)
    ii <- ii+1
    if (ppc.options$debug) cat(ii)
    pk1[[ii]] <- cbind(logmz[aa],a[,2])
  }
  
  return(pk1)
}
