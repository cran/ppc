ppc.error<- function(y,yhat){
  ## compute errors and confusion matrices, from a ppc fit
  
  tt <- vector("list",ncol(yhat))
  codesy<-levels(y)
  yhat2<-matrix(codesy[yhat], nrow=nrow(yhat),ncol=ncol(yhat))
  
  
  for(i in 1:length(tt)){
    tt[[i]] <- table(y,yhat2[,i])
    
  }
  err <- rep(NA,length(tt))
  
  for(i in 1:length(tt)){
    err[i] <- sum(yhat2[,i]!=y)
  }
  
  return(list(err=err,confusion=tt))
}


