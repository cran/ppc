ppc.cv <- function(ppc.fit, data,   user.parms){
  ## K-fold cross-validation for ppc method
  ## takes results of make.centroids.list, find.splits, and predict.ppc
  ## and computes cross-validated predictions and error estimate.
  
  peaklist.tr<- data$peaklist
  ytr <- data$ytr
  logmz <- data$logmz
  
  peak.gap <- user.parms$peak.gap
  nsplits <- user.parms$nsplits
  fix.at.one <- user.parms$fix.at.one
  recluster <- user.parms$recluster
  
  folds <- balanced.folds(ytr)
  
  n.class<-table(ytr)
  n.threshold<- length(ppc.fit$threshold)
  
  n <- length(peaklist.tr)
  yhatcv <- array(NA,c(n,n.threshold))
  probcv <- array(1, c(n, length(n.class), n.threshold))
  
  for(ii in 1:length(folds)){
    gg <- folds[[ii]]
    if (ppc.options$debug) cat(c("fold=",ii),fill=T)
    
    pk <- peaklist.tr[-gg]
    
    data.temp <- list(ytr=ytr[-gg], logmz=logmz, peaklist=pk)
    astar <- ppc.make.centroid.list(data.temp,  user.parms)
    
    junk0 <- ppc.predict.peaks(astar,  data.temp)
    aa <- ppc.find.splits(astar, junk0 ,data.temp,user.parms)
    ress <- ppc.predict(astar,aa,logmz,peaklist.te=peaklist.tr[gg],threshold=ppc.fit$threshold, metric="euclidean")
    yhatcv[gg,] <- ress$yhat
    probcv[gg,,] <- ress$prob
  }
  
  junk <- ppc.cv.error(ytr, yhatcv, folds)
  
  return(list(err=junk$err,
              se=junk$se,
              confusion=junk$confusion,
              threshold=ppc.fit$threshold,
              yhat=yhatcv,
              prob=probcv,
              y=ytr,
              folds=folds,
              numsites=ppc.fit$numsites))
}
