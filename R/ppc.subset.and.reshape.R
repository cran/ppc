##
## Original version
##

#ppc.subset.and.reshape<- function(data, user.parms){
  
#  xtr<-data$xtr
#  mz<-data$mz
#  logmz<-data$logmz
#  ytr<-data$ytr
  
  
#  mz.labels<-NULL
#  sample.labels <- data$sample.labels
  
#  if (!is.null(user.parms$mz.min) | !is.null(user.parms$mz.max)) {
#    o1<- user.parms$mz.min
#    o2<- user.parms$mz.max
#    if(!is.null(o1)) {
#      xtr<- xtr[mz>= o1,]
#      logmz<- logmz[mz>= o1]
#      mz<- mz[mz>= o1]
#    }
#    if (!is.null(o2)) {
#      xtr<- xtr[mz<= o2,]
#      logmz<- logmz[mz<= o2]
#      mz<-    mz[mz<= o2]
#    }
#  }
#  ## here I reshape the data into one column per patient. mz and logmz values are
#  ## string out into one long vector, spaced apart

#  if (!is.null(data$batch.labels)) {
#    npatients<-length(unique(data$patient.labels))
#    nbatches<-length(unique(data$batch.labels))
#    p<-nrow(xtr)
    
#    labelord<-order(data$batch.labels)
#    xtr<-xtr[,labelord]
#    if(!is.null(ytr)){ ytr<- ytr[labelord][1:npatients]}
    
#    xtr3<-matrix(NA,nrow=nbatches*nrow(xtr),ncol=npatients)
#    ii<-1
#    for (i in 1:nbatches) {
#      xtr3[ii:(ii+p-1),]<- xtr[, ((i-1)*npatients+1):(i*npatients)]
#      ii<-(ii+p)
#    }
    
    
#    mz4<-NULL
#    logmz<-NULL
#    mz.labels<-NULL
#    fac<- 2*(max(mz)-min(mz))
#    faclog<- 2*(max(log(mz))-min(log(mz)))
    
#    for(ii in 1:nbatches) {
#      mz4<-c(mz4,mz+(ii-1)*fac)
#      logmz<-c(logmz,log(mz)+(ii-1)*faclog)
#      mz.labels<-c(mz.labels,paste(round(mz,1),"-",batches[ii],sep=""))
#    }
#    xtr<-xtr3
#    mz<-mz4
    
#    sample.labels<-data$patient.labels[labelord]
#    sample.labels<-sample.labels[1:npatients]
#  }
  
#  iytr<-NULL
#  if(!is.null(ytr)) {
#    iytr<- as.factor(ytr)
#  }
  
#  return(list(mz=mz,logmz=logmz, mz.labels=mz.labels, 
#              xtr=xtr, ytr=iytr, sample.labels=data$sample.labels))
#}

ppc.subset.and.reshape<- function(data, user.parms){
  
 data.keep<-data

  xtr<-data$xtr
  mz<-data$mz
  logmz<-data$logmz
  ytr<-data$ytr
  
 mz.labels<-NULL

  sample.labels <- data$sample.labels
  patient.labels<- data$patient.labels
 
  if (!is.null(user.parms$mz.min) | !is.null(user.parms$mz.max)) {
    o1<- user.parms$mz.min
    o2<- user.parms$mz.max
    if(!is.null(o1)) {
      xtr<- xtr[mz>= o1,]
      logmz<- logmz[mz>= o1]
      mz<- mz[mz>= o1]
    }
    if (!is.null(o2)) {
      xtr<- xtr[mz<= o2,]
      logmz<- logmz[mz<= o2]
      mz<-    mz[mz<= o2]
    }
  }
  ## here I reshape the data into one column per patient. mz and logmz values are
  ## string out into one long vector, spaced apart

  if (!is.null(data$batch.labels)) {
    npatients<-length(unique(data$patient.labels))
    nbatches<-length(unique(data$batch.labels))
    p<-nrow(xtr)
    
    labelord<-order(data$batch.labels)
    xtr<-xtr[,labelord]
    if(!is.null(ytr)){ ytr<- ytr[labelord][1:npatients]}
    
    xtr3<-matrix(NA,nrow=nbatches*nrow(xtr),ncol=npatients)
    ii<-1
    for (i in 1:nbatches) {
      xtr3[ii:(ii+p-1),]<- xtr[, ((i-1)*npatients+1):(i*npatients)]
      ii<-(ii+p)
    }
    
    
    mz4<-NULL
    logmz<-NULL
# ROB fixed this bug (subtle error for pred test data)

    mz.labels<-NULL
#    fac<- 2*(max(mz)-min(mz))
#    faclog<- 2*(max(log(mz))-min(log(mz)))
    
fac<- 2*(user.parms$mz.max- user.parms$mz.min)
faclog<- 2*(log(user.parms$mz.max)- log(user.parms$mz.min+1))

    for(ii in 1:nbatches) {
      mz4<-c(mz4,mz+(ii-1)*fac)
      logmz<-c(logmz,log(mz)+(ii-1)*faclog)
      mz.labels<-c(mz.labels,paste(round(mz,1),"-",data$batches[ii],sep=""))
    }
    xtr<-xtr3
    mz<-mz4
    
    sample.labels<-data$patient.labels[labelord]
    sample.labels<-sample.labels[1:npatients]
    patient.labels<-sample.labels
   
  }
  
  iytr<-NULL
  if(!is.null(ytr)) {
    iytr<- as.factor(ytr)
  }
  
 data.keep$mz<-mz
 data.keep$logmz<-logmz
 data.keep$mz.labels<- mz.labels
 data.keep$xtr<-xtr
 data.keep$ytr<-iytr
 data.keep$sample.labels <- sample.labels
 data.keep$patient.labels <- patient.labels
 
 
  return( data.keep)
}


