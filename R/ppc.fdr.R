##
## Change by naras. Added split.fit argument.
##
ppc.fdr <- function(data, centroid.fit, peak.fit, split.fit, ppc.fit, user.parms){
  ytr <- data$ytr
  fix.at.one <- user.parms$fix.at.one
  nsplits <- user.parms$nsplits
  nperms <- user.parms$nperms
  m <- nrow(centroid.fit$cen)
  
  threshold <- ppc.fit$threshold
  tt <- split.fit$prhat[,1]-split.fit$prhat[,2]
  
  ress <- vector("list",nperms)
  ttstar <- matrix(NA,nrow=length(tt),ncol=nperms)
  
  for(i in 1:nperms){
    
    if (ppc.options$debug) cat(c("i=",i),fill=T)
    ytr2 <- sample(ytr)
    data2 <- data
    data2$ytr <- ytr2
    split.fit2 <- ppc.find.splits(centroid.fit, peak.fit, data2, user.parms) 
    
    ress[[i]] <- split.fit2$prhat
    ttstar[,i] <- split.fit2$prhat[,1]-split.fit2$prhat[,2]
  }
  
  nt <- length(threshold)
  npeaks <- rep(NA,nt)
  npeaks0 <- npeaks
  for(i in 1:nt){
    npeaks[i] <- sum(abs(split.fit$prhat[,1]-split.fit$prhat[,2])>threshold[i])
    npeaks0[i] <- 0
    for(j in 1:nperms){
      npeaks0[i] <- npeaks0[i]+sum(abs(ress[[j]][,1]-ress[[j]][,2])>threshold[i])
    }}
  
  fdr <- (npeaks0/nperms)/npeaks
  q1 <- quantile(abs(tt), .25)
  q2 <- quantile(abs(tt), .75)
  
  pi0 <- min((sum(abs(ttstar)> q1 & abs(ttstar)< q2)/nperms)/(.5*m) ,1 )
  
  fdr <- fdr*pi0
  
  results <- cbind(threshold, npeaks,fdr)
  
  dimnames(results) <- list(NULL,c("threshold", "npeaks","fdr"))
  return(list(results=results,pi0=pi0, threshold=threshold))
}
