##
## Original by Rob
##
#ppc.subset.and.reshape.peakdata<- function(data, user.parms){


## subset and reshaping for data containing initial peak lists

                      
#peaklist<-data$peaklist
#logmz<-data$logmz
#mz<-data$mz
#ytr<-data$ytr

#mz.labels<-NULL
#n<-length(peaklist)

#if(!is.null(user.parms$mz.min) | !is.null(user.parms$mz.max)){
#	 o1<- user.parms$mz.min
#	 o2<- user.parms$mz.max
#	 if(!is.null(o1)){ 
#                      peaklist<-subset.peaklist(peaklist,log(o1),-1)
#                   logmz<- logmz[mz>= o1]
#                   mz<- mz[mz>= o1]
                  
#  	}
#  	if(!is.null(o2)){
#                      peaklist<-subset.peaklist(peaklist,log(o2),+1)
#                     logmz<- logmz[mz<= o2]
#                     mz<-    mz[mz<= o2]
#  	}
#}

## here I reshape the data into peaklist per patient. mz and logmz values are
## strong out into one long vector, spaced apart

#if(!is.null(data$batch.labels)){
#         nbatches<- length(unique(data$batch.labels))
#        npatients<-n/nbatches

#	p<-length(peaklist)

#	patient.ord<-order(data$patient.labels)
       
#	pk<-peaklist[patient.ord]
      
#         mz4<-NULL
#        logmz<-NULL
#        mz.labels<-NULL
#        fac<- 2*(max(mz)-min(mz))
#        faclog<- 2*(max(log(mz))-min(log(mz)))

#       pknew<-vector("list",npatients)
#        ii<-0
#	for(i in 1:npatients){
#           temp0<-NULL
#           for(j in 1:nbatches){
#             ii<-ii+1
#             temp<-pk[[ii]]
#             temp[,1]<-temp[,1]+faclog*(j-1)
#	     temp0<-rbind(temp0,temp)
#           }
#           pknew[[i]]<-temp0
#	}
	

#	mz.labels<-NULL
#        fac<- 2*(max(mz)-min(mz))
#        faclog<- 2*(max(log(mz))-min(log(mz)))

#	for(ii in 1:nbatches){
#	   mz4<-c(mz4,mz+(ii-1)*fac)
#	  logmz<-c(logmz,log(mz)+(ii-1)*faclog)
#	  mz.labels<-c(mz.labels,paste(round(mz,1),"-",data$batches[ii],sep=""))
#        }
#mz<-mz4
#if(!is.null(ytr)){  ytr<- matrix(ytr[patient.ord],ncol=nbatches,byrow=T)[,1]}
#sample.labels<- matrix(data$sample.labels[patient.ord],ncol=nbatches,byrow=T)[,1]
#peaklist<- pknew
#}

#iytr<-NULL
#if(!is.null(ytr)){ iytr<- as.factor(ytr)}

#return(list(mz=mz,logmz=logmz, mz.labels=mz.labels, peaklist=peaklist, ytr=iytr,
#sample.labels=sample.labels))
#}


"subset.peaklist"<- function(peaklist,r, dir){
       n<-length(peaklist)
       for(i in 1:n){
        if(dir== -1){
           o<-peaklist[[i]][,1]>r
              peaklist[[i]]<-peaklist[[i]][o,]
        }
     if(dir== +1){
           o<-peaklist[[i]][,1]<r
              peaklist[[i]]<-peaklist[[i]][o,]
        }
}
return(peaklist)
}
 ppc.subset.and.reshape.peakdata<-
function (data, user.parms) 
{
  data.keep <- data
    peaklist <- data$peaklist
    logmz <- data$logmz
    mz <- data$mz
    ytr <- data$ytr
    sample.labels <- data$sample.labels
    mz.labels <- NULL
    n <- length(peaklist)
    if (!is.null(user.parms$mz.min) | !is.null(user.parms$mz.max)) {
        o1 <- user.parms$mz.min
        o2 <- user.parms$mz.max
        if (!is.null(o1)) {
            peaklist <- subset.peaklist(peaklist, log(o1), -1)
            logmz <- logmz[mz >= o1]
            mz <- mz[mz >= o1]
        }
        if (!is.null(o2)) {
            peaklist <- subset.peaklist(peaklist, log(o2), +1)
            logmz <- logmz[mz <= o2]
            mz <- mz[mz <= o2]
        }
    }
    if (!is.null(data$batch.labels)) {
        nbatches <- length(data$batches)
        npatients <- n/nbatches
        p <- length(peaklist)
        patient.ord <- order(data$patient.labels)
        pk <- peaklist[patient.ord]
        mz4 <- NULL
        logmz <- NULL
        mz.labels <- NULL

# ROB fixed this bug (subtle error for test set pred!)

#        fac <- 2 * (max(mz) - min(mz))
#        faclog <- 2 * (max(log(mz)) - min(log(mz)))

        fac<- 2*(user.parms$mz.max- user.parms$mz.min)
        faclog<- 2*(log(user.parms$mz.max)- log(user.parms$mz.min+1))

        pknew <- vector("list", npatients)
        ii <- 0
        for (i in 1:npatients) {
            temp0 <- NULL
            for (j in 1:nbatches) {
                ii <- ii + 1
                temp <- pk[[ii]]
                temp[, 1] <- temp[, 1] + faclog * (j - 1)
                temp0 <- rbind(temp0, temp)
            }
            pknew[[i]] <- temp0
        }
        mz.labels <- NULL
        fac <- 2 * (max(mz) - min(mz))
        faclog <- 2 * (max(log(mz)) - min(log(mz)))
        for (ii in 1:nbatches) {
            mz4 <- c(mz4, mz + (ii - 1) * fac)
            logmz <- c(logmz, log(mz) + (ii - 1) * faclog)
            mz.labels <- c(mz.labels, paste(round(mz, 1), "-", 
                data$batches[ii], sep = ""))
        }
        mz <- mz4
        if (!is.null(ytr)) {
            ytr <- matrix(ytr[patient.ord], ncol = nbatches, 
                byrow = T)[, 1]
        }
        sample.labels <- matrix(data$sample.labels[patient.ord], 
            ncol = nbatches, byrow = T)[, 1]
        peaklist <- pknew
    }
    iytr <- NULL
    if (!is.null(ytr)) {
        iytr <- as.factor(ytr)
    }
    data.keep$mz <- mz
    data.keep$logmz <- logmz
    data.keep$mz.labels <- mz.labels
    data.keep$peaklist <- peaklist
    data.keep$ytr <- iytr
    data.keep$sample.labels <- sample.labels
    return(data.keep)
}

