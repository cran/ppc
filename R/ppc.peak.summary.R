ppc.peak.summary <- function(centroid.fit, peak.fit, data, user.parms, split.fit) {
  temp<-  peak.fit$ind*peak.fit$ht
  temp[temp==0] <- NA
  
  sum.na<- function(x){sum(x[!is.na(x)])}
  
  splits <- rep(NA,length(split.fit$cuthat))
  for(i in 1:length(split.fit$cuthat)){
    splits[i]<- split.fit$cutpoints[i,split.fit$cuthat[i]]
  }
 nc<-ncol(split.fit$prhat)

  labs0 <- rep(NA,nc)
  for(i in 1:nc){
    labs0[i] <- paste("Prop.in.class.",as.character(i),sep="")
}
  prdiff<-matrix(NA,ncol=nc*(nc-1)/2,nrow=nrow(split.fit$prhat))
  ii<-0
  labs <- rep(NA,nc*(nc-1)/2)
  for(i in 1:(nc-1)){
    for(j in (i+1):nc){
      ii <- ii+1
  prdiff[,ii] <- split.fit$prhat[,j]-split.fit$prhat[,i]
      labs[ii] <- paste("Pr",as.character(j),"-","Pr",as.character(i),sep="")
    }}
      
  rank.split <- rank(-apply(abs(prdiff),1,sum))
  
  res <- cbind(exp(centroid.fit$cent[,-2]), apply(temp>0,1,sum.na), rank.split,splits,split.fit$prhat, prdiff, temp)
  
  dimnames(res) <- list(NULL,c("peak.position",
                               "min",
                               "max",
                               "number.of.spectra",
                               "rank.of.peak",
                               "height.split.point",
                               labs0,
                               labs,
                               data$sample.labels))
  return(res)
}



