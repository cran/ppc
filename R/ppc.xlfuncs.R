##
## Functions for the Excel interface
##
##
## Read all data files in a folder and return mz and the matrix of spectra
##

##
## A list that is used for setting ppc.options
##
ppc.options <- list(debug=FALSE, #whether to turn on debugging or not
                    err.file=ifelse(.Platform$OS.type=="windows", "C:/ppctrace.txt", "ppctrace.txt"),
                    reserved.class.label="Unspecified")

ppc.constants <- list(data.types = list(raw="RAW", peak="PEAK"))


##
## Our error handler
##
ppc.xl.error.trace <- function() {
  err.message <- geterrmessage()
  sink(ppc.options$err.file)
  print(err.message)
  traceback()
  sink()
  winDialog(type="ok", message=err.message)
}

##
## Upon loading, if we are in a windows environment, we use the windows
## dialog mechanism to display errors. Useful for debugging COM apps
##
.onLoad <- function(lib, pkg) {
  if ( .Platform$OS.type == "windows") {
#    options(error=function() winDialog(type="ok", message=geterrmessage()))
    options(error=ppc.xl.error.trace)
  }

}

##
## Upon unload, we set things back the way they were...
##
.onUnload <- function(libpath){
  if ( .Platform$OS.type == "windows") {
    options(error=NULL)
  }
}
       
##
## A function to get the names of subdirs in a given directory.
## Note: full pathnames are returned.
##
ppc.xl.get.names.of.subdirs <- function(directory) {
  f.info <- file.info(dir(directory, full.names=TRUE))
  return(dimnames(f.info)[[1]][f.info[, "isdir"] == TRUE])
}

##
## A function to get the names of files in a given directory.
## Note: full pathnames are returned.
##
ppc.xl.get.names.of.files <- function(directory) {
  f.info <- file.info(dir(directory, full.names=TRUE))
  return(dimnames(f.info)[[1]][f.info[, "isdir"] == FALSE])
}

##
## NOTE to memory-leaking myself:
##       All the functions below expect the class data directory
##       that is, the directory that is the root of all class specific data
##       We have to call the functions for each class. The plurals expect the entire list.
##

##
## A function to get detect if we have batches or not.
##
ppc.xl.detect.batches <- function(classDataDirectory) {
  x <- ppc.xl.get.names.of.subdirs(classDataDirectory)
  return(length(x) >= 2)  # Need at least two batches!
}

##
## A function to get detect the batch labels
## Note: only the batch labels (not full pathname) are returned
##       Also, batches are assumed to be present.
##       Batch labels are merely base names of subdirectories 
##
ppc.xl.compute.batch.labels <- function(classDataDirectories) {
  x <- unlist(sapply(classDataDirectories,
                     function(x) ppc.xl.get.names.of.files(ppc.xl.get.names.of.subdirs(x))))
  return (as.vector(sapply(x, function(x) basename(dirname(x)))))
}

##
## A function to get the list of all data files 
## Note: Full path names are returned
##
ppc.xl.get.all.data.file.names <- function(classDataDirectories, batches.present=FALSE) {
  if (batches.present) {
    x <- unlist(sapply(classDataDirectories,
                       function(x) ppc.xl.get.names.of.files(ppc.xl.get.names.of.subdirs(x))))
  } else {
    x <- unlist(sapply(classDataDirectories,
                       function(x) ppc.xl.get.names.of.files(x)))
  }
  return(as.vector(x))
}

##
## A function to compute sample labels
## Note: if batches are present, sample labels are <batch_name>-<file-name-without-suffix>
##       otherwise just <file-name-without-suffix>
##
ppc.xl.compute.sample.labels <- function(classDataDirectories, batches.present=FALSE) {
  if (batches.present) {
    x <- unlist(sapply(classDataDirectories,
                       function(x)
                       basename(ppc.xl.get.names.of.files(ppc.xl.get.names.of.subdirs(x)))))
    x <- paste(ppc.xl.compute.batch.labels(classDataDirectories), x, sep="-")
  } else {
    x <- unlist(sapply(classDataDirectories,
                       function(x)
                       basename(ppc.xl.get.names.of.files(x))))
  }
  return(as.vector(sapply(x, ppc.xl.remove.suffix)))
}

##
## A function to compute sample sizes
## Note: Input is the vector of data folders
##
ppc.xl.compute.sample.sizes <- function(classDataDirectories, batches.present=FALSE) {
  result <- sapply(classDataDirectories,
                   function(x)
                   ifelse(batches.present,
                          length(ppc.xl.get.names.of.files(ppc.xl.get.names.of.subdirs(x))),
                          length(ppc.xl.get.names.of.files(x))))
  return (as.vector(result))
}


##
## A function to compute patient labels
## Note: patient labels are merely basenames of files without the extension
##
ppc.xl.compute.patient.labels <- function(classDataDirectories, batches.present=FALSE) {
  if (batches.present) {
    x <- unlist(sapply(classDataDirectories,
                       function(x)
                       basename(ppc.xl.get.names.of.files(ppc.xl.get.names.of.subdirs(x)))))
  } else {
    x <- unlist(sapply(classDataDirectories,
                       function(x)
                       basename(ppc.xl.get.names.of.files(x))))
  }
  return(as.vector(sapply(x, ppc.remove.suffix)))
}


##
## A function to get the names of files in a given directory for a given batch
## Note: full pathnames are returned
##
ppc.xl.get.files.in.batch <- function(classDataDirectory, batchName) {
  f.info <- file.info(dir(file.path(classDataDirectory, batchName), full.names=TRUE))
  return(dimnames(f.info)[[1]][f.info[, "isdir"] == FALSE])
}



##
## Read data from a single file. 
##
## 
ppc.xl.read.data.file  <- function(file.name, mz=NULL) {
  if (ppc.options$debug) {
    print(paste("Reading file", file.name))
  }

  file.contents <- read.table(file.name, sep=",")

  if (!is.null(mz)) {
    y <- approx(file.contents[, 1], file.contents[, 2], xout = mz)$y
  } else {
    mz <- file.contents[, 1]
    y <- file.contents[, 2]
  }
  return(list(mz=mz, y = y))
}

##
## Build raw data set
##
ppc.xl.build.data  <- function(raw.file.data,
                               class.labels,
                               sample.labels,
                               batch.labels,
                               patient.labels,
                               batches.exist = FALSE,
                               data.type=ppc.constants$data.types$raw,
                               class.levels = NULL) {
  if (is.null(class.levels)) {
    ytr <- factor(class.labels)
  } else {
    ytr <- factor(class.labels, levels=class.levels)
  }

  batches <- NULL
  if (batches.exist) {
    batches <- unique(batch.labels)
  }
  
  if (data.type == ppc.constants$data.types$raw) {
    return(list(mz=raw.file.data$mz,
                logmz = log(raw.file.data$mz),
                xtr = raw.file.data$xtr,
                ytr = ytr,
                peaklist = raw.file.data$peaklist, ## will be NULL!
                sample.labels = sample.labels,
                batch.labels = batch.labels,
                batches = batches,
                patient.labels = patient.labels,
                data.type = data.type))
  } else {
    return(list(mz=raw.file.data$mz,
                logmz = log(raw.file.data$mz),
                xtr = NULL,
                ytr = ytr,
                peaklist = raw.file.data$peaklist,
                sample.labels = sample.labels,
                batch.labels = batch.labels,
                batches = batches,
                patient.labels = patient.labels,
                data.type = data.type))
  }
}

##
## Get the current user parameters
##
ppc.xl.get.user.parameters  <- function() {
  return(ppc.xl.current.user.parameters)
}

##
## Set the current user parameters
##
ppc.xl.set.user.parameters  <- function(x) {
  ppc.xl.current.user.parameters  <- x
}

##
## Return the default set of parameters
##
ppc.xl.get.default.parameters  <- function() {
  return (list(stn = 1,
               minht = 0,
               span = 201,
               smoothing.span = 0.05,
               nsplits = 10,
               fix.at.one = FALSE,
               peak.gap = 0.005,
               recluster = FALSE,
               mz.min = 15,
               mz.max = 1500,
               nperms=20))
}

##
## Construct the current set of parameters
##
ppc.xl.make.parameter.set  <- function(stn=1, minht=0, span=201,
                                       smoothing.span=0.05,
                                       nsplits=10,
                                       fix.at.one=FALSE,
                                       peak.gap=0.005,
                                       recluster=FALSE,
                                       mz.min=15,
                                       mz.max=1500,
                                       nperms=20) {
  return (list(stn = stn,
               minht = minht,
               span = span,
               smoothing.span = smoothing.span,
               nsplits = nsplits,
               fix.at.one = fix.at.one,
               peak.gap = peak.gap,
               recluster = recluster,
               mz.min = mz.min,
               mz.max = mz.max,
               nperms=nperms))
}


##
## Compute the training errors for the dataset
ppc.xl.compute.training.errors  <- function(fit, data) {
  foo  <- ppc.error(data$ytr, fit$yhat)
  threshold = fit$threshold
  n  <- length(data$ytr)
  
  return (list(x = fit$threshold,
               y = foo$err/n,
               y.ytop = fit$numsites,
               x.label = "Threshold",
               y.label = "Training Error"))
}

##
## Compute the test errors for the dataset
## TODO Need to fix.
ppc.xl.compute.test.errors  <- function(fit, data) {
  err.cnts  <- ppc.error(data$ytr, fit$yhat)
  threshold <- fit$threshold
  n  <- length(data$ytr)
  
  return (list(x = fit$threshold,
               y = err.cnts$err/n,
               y.ytop = fit$numsites,
               x.label = "Threshold",
               y.label = "Test Error"))
}


##
## Rob likes to predict a whole set of things in one shot, but there is nothing that says
## that a user might not want to predict for a particular value of a threshold.
## SO I have to rewrite the predict.ppc function to do things just for a single threshold.
##
ppc.xl.predict <- function(centroids.fit, split.fit, logmz, peaklist.te, threshold=NULL,
                           metric=c("binomial","euclidean", "absolute"),
                           summ=c("mean","median")) {
  
  ## test set prediction for PPC method
  ## makes  predictions  for test set list of peaks in peaklist.te
  ##  
  ## "metric" is the metric used in the nearest centroid rule; "summ" is
  ## summary used  to combine distances over sites
  ## n.threshold is the number of shrinkage thresholds used
  
  ## returns yhat, plus ind0= indicator matrix  of whether peak was found,
  ##           ind=indicator matrix  of whether peak > cutpoint was found
  ##           ht=matrix of peak heights
  
  ## "numsites" is the number of sites present at each threshold; "sites" is a list
  ## containing the ##s of these sites
  ##
  ## dis is the distance of each test profile to each class centroid
  
  ##
  metric <- match.arg(metric)
  summ <- match.arg(summ)
  ##  print(metric)
  ##  print(summ)
  
  n <- length(peaklist.te)
  n.class <- split.fit$n.class
  K <- ncol(split.fit$prhat)
  p <-length(split.fit$cuthat)
  yhat <- matrix(NA, nrow = p, ncol = n)
  yhatt <- vector(mode="numeric", length = n.class)
  ind0 <- yhat
  ht <- yhat
  dis <- array(NA,c(n, K))

  for(j in 1:n) {
    if (ppc.options$debug) cat(j)
    aaa <- ppc.predict.peaks1(centroids.fit, logmz,  peaklist= 
                              peaklist.te[[j]])
    ht.hat <-  aaa$ht * aaa$ind
    ind0[,j] <- 1*aaa$ind
    ht[,j] <- ht.hat
    for(i in 1:p) {
      yhat[i, j] <- 1 * (ht.hat[i] > split.fit$cutpoints[i, split.fit$cuthat[i]])
    }
  }        
  
  prmean <- apply(split.fit$prhat,1,mean)
  delta <- split.fit$prhat-matrix(prmean,nrow=nrow(split.fit$prhat),ncol=K)
  
  dd<- matrix(NA,nrow=n,ncol=K)
  
  ##  if(is.null(threshold)){threshold <- seq(0, max(abs(delta)), length = n.threshold)}
  
  soft.thresh <- function(x, tt) {
    sign(x) * (abs(x) - tt) * (abs(x) > tt)
  }
  
  delta2 <- soft.thresh(delta, threshold)
  pr <- prmean + delta2
  sumabs<-  apply(abs(delta2),1,sum)
  pos <-  sumabs!=0
  numsites <- sum(pos)
  sites <- (1:nrow(pr))[pos]
  
  sites <- sites[order( - abs(sumabs[pos]))]
  if(metric == "euclidean") {
    for(k in 1:K){
      temp<- (yhat - matrix(pr[,k], ncol = n, nrow = nrow(yhat)))^2
      dd[,k] <- apply(temp,2,summ)
    }
  }
  if(metric == "absolute") {
    for(k in 1:K){
      temp<- abs(yhat - matrix(pr[,k], ncol = n, nrow = nrow(yhat)))
      dd[,k] <- apply(temp,2,summ)
    }
  }
  if(metric == "binomial") {
    ## compute Bayes estimate of pr, to keep it away from 0 or 1
    ntemp <- matrix(n.class,nrow=nrow(pr),ncol=length(n.class), byrow=T)
    prbayes<- (pr*ntemp+1)/(ntemp+2)
    for(k in 1:K){
      prmat<-matrix(prbayes[,k], ncol = n, nrow = nrow(yhat))
      temp<- -1*(yhat*log(prmat)+(1-yhat)*log(1-prmat))
      dd[,k] <- apply(temp,2,summ)
    }
  }
  
  yhatt <- apply(dd,1,which.is.min)
  dis <- dd

  prob <- exp(-dis)
  prob <- prob/apply(prob,1,sum)
  
  return(list(yhat = yhatt, threshold = threshold, numsites = numsites,
              sites = sites, ind = yhat, ind0=ind0, prob=prob, ht=ht))
  
}


##
## Compute the confusion matrix for the dataset
ppc.xl.compute.confusion.matrix  <- function(centroids.fit, split.fit,
                                             data, peaklist.te, threshold=NULL,
                                             metric=c("binomial","euclidean", "absolute"),
                                             summ=c("mean","median")) {
  
  ##i has to be determined!
  metric  <- match.arg(metric)
  summ  <- match.arg(summ)
  pred.obj  <- ppc.xl.predict(centroids.fit, split.fit,
                              data$logmz, peaklist.te=peaklist.te, threshold=threshold,
                              metric=metric, summ=summ)
  
  mat  <- table(data$ytr, pred.obj$yhat)
  return(list(confusion.matrix=mat))
}

##
## massage arrays and matrices to handle NAs in Excel
##
ppc.xl.transform.matrix  <- function(x) {
  w  <- x
  m  <- apply(w, 2, function(x) {ifelse(is.na(x),1,0)})
  w  <- apply(w, 2, function(x) {ifelse(is.na(x),0,x)})
  return(list(matrix=w, missing=m))
}

##
## Compute the ppc cv errors
##
ppc.xl.compute.cv.errors  <- function(trained.obj, cv.obj) {
  n <- nrow(cv.obj$yhat)
  y <- cv.obj$y
  nc <- length(table(y))
  nfolds <- length(cv.obj$folds)
  err2 <- matrix(NA, nrow = length(unique(y)), ncol = length(cv.obj$threshold))
  for(i in 1:(length(cv.obj$threshold))) {
    s <- cv.obj$confusion[[i]]
    diag(s) <- 0
    err2[, i] <- apply(s, 1, sum)/table(y)
  }
  
  return(list(x=trained.obj$threshold, y=cv.obj$err,
              x.label="Threshold", y.label="Error",
              y.se=cv.obj$se, numsites=cv.obj$numsites,
              cv.err=t(err2), cv.legend= dimnames(table(y))[[1]]))
}

##
## Find the nearest value to a list of nearest value
##
ppc.xl.find.index.of.nearest.value  <- function(target.vec, val) {
  return(which(order(abs(target.vec - val)) == 1))
}


##
## Setup for cv
##
ppc.xl.cv.setup <- function(ppc.fit, data,   user.parms) {
  
  ## K-fold cross-validation for ppc method
  ## takes results of make.centroids.list, find.splits, and predict.ppc
  ## and computes cross-validated predictions and error estimate.

  peaklist.tr <- data$peaklist
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
  return(list(numsites = ppc.fit$numsites,
              peaklist.tr = peaklist.tr,
              ytr = ytr,
              logmz = logmz,
              peak.gap = peak.gap,
              nsplits = nsplits,
              fix.at.one = fix.at.one,
              recluster = recluster,
              folds = folds,
              threshold = ppc.fit$threshold,
              n.threshold = n.threshold,
              n = n,
              yhatcv = yhatcv,
              probcv = probcv))
}

##
## DO fold ii
##
ppc.xl.cv.do.fold <- function(cv.state, ii, user.parms) {
  gg <- cv.state$folds[[ii]]
  if (ppc.options$debug) print(c("fold=",ii),fill=T)
  pk <- cv.state$peaklist.tr[-gg]
  data.temp <- list(ytr=cv.state$ytr[-gg], logmz=cv.state$logmz, peaklist=pk)
  astar <- ppc.make.centroid.list(data.temp,  user.parms)
  junk0 <- ppc.predict.peaks(astar,  data.temp)
  aa <- ppc.find.splits(astar, junk0 ,data.temp,user.parms)
  ress <- ppc.predict(astar,aa,cv.state$logmz,peaklist.te=cv.state$peaklist.tr[gg],
                      threshold=cv.state$threshold, metric="euclidean")
  return(list(gg=gg, yhat=ress$yhat, prob=ress$prob))
  ##  yhatcv[gg,] <- ress$yhat
  ##  cv.state$probcv[gg,,] <- ress$prob
}

##
## Compute the cv probabilities
##
ppc.xl.compute.cvprobs  <- function(cv.fit, data, threshold) {
  ii <- (1:length(cv.fit$threshold))[cv.fit$threshold > threshold]
  ii <- ii[1]
  ss <- data$sample.labels
  pp <- cv.fit$prob[,  , ii]
  
  y <- cv.fit$y
  
  o <- order(y)
  y <- y[o]
  if(!is.null(ss)) {
    ss <- ss[o]
  }
  ppp <- pp[o,  ]
  n <- nrow(ppp)
  nc <- length(unique(y))
  
  return (list(x = 1:n,
               y = ppp,
               x.label = "Sample",
               y.label = "CV Probabilities",
               y.names = levels(factor(y)),
               y.lines = cumsum(table(cv.fit$y)),
               x.dummy = vector(length=2, mode="numeric"),
               y.dummy = vector(length=2, mode="numeric"),
               x.names = ss))
}


##
## Compute the prediction probabilities for test set.
##
  
ppc.xl.compute.test.probs  <- function(centroid.fit, split.fit, training.class.names, data, peaklist,
                                       threshold) {
  ppc.fit1 <-  ppc.predict1(centroid.fit, split.fit, data$logmz, peaklist, threshold)
  predicted.y <- sapply(ppc.fit1$yhat, function(x){ training.class.names[x] })
  predicted.probs <- ppc.fit1$prob
  sample.labels <- data$sample.labels

  order.classes  <- order(predicted.y)
  actual.classes <- as.character(data$ytr[order.classes])
  actual.classes[is.na(actual.classes)] <- ppc.options$reserved.class.label
  pp <- apply(predicted.probs, 2, function(x) x[order.classes])
  ny  <- predicted.y[order.classes]
  n  <- length(ny)
  sample.labels  <- data$sample.labels[order.classes]
  actual.class.names <- levels(factor(actual.classes))

  return (list(x = 1:n,
               y = pp,
               x.label = "Sample",
               y.label = "Predicted Test Probabilities",
               y.names = training.class.names,
               y.lines = cumsum(table(actual.classes)) + 0.5,
               x.dummy = vector(length=2, mode="numeric"),
               y.dummy = vector(length=2, mode="numeric"),
               panel.names = actual.class.names,
               x.names = sample.labels))
}

##
## Compute the quantities for the ``Show Prediction'' button in PPC (Excel)
##

ppc.xl.compute.test.prediction  <- function(centroid.fit, split.fit, training.class.names, data, peaklist,
                                       threshold) {
  ppc.fit1 <-  ppc.predict1(centroid.fit, split.fit, data$logmz, peaklist, threshold)
  predicted.y <- sapply(ppc.fit1$yhat, function(x){ training.class.names[x] })
  predicted.probs <- ppc.fit1$prob
  colnames(predicted.probs) <- training.class.names
  actual.classes <- as.character(data$ytr)
  data.has.missing.class.labels <- any(is.na(actual.classes))
  actual.classes[is.na(actual.classes)] <- ppc.options$reserved.class.label
  
  if (data.has.missing.class.labels) {
    confusion.matrix <- NULL
  } else {
    confusion.matrix <- table(actual.classes, predicted.y)
  }
  return (list(actual.class.labels = actual.classes,
               predicted.class.labels = predicted.y,
               sample.labels = data$sample.labels,
               patient.labels = data$patient.labels,
               confusion.matrix = confusion.matrix,
               predicted.probs = t(predicted.probs)))
}



##
## Transform class labels by replacing NAs with reserved label
##
ppc.xl.transform.class.labels <- function(class.labels, reservedLabel) {
  w <- class.labels
  w[is.na(class.labels)] <- reservedLabel
  return (w)
}

##
## Plotting of histograms
##
ppc.xl.plothist <-  function(peak.fit, ppc.fit, centroid.fit, split.fit, data, nsites=NULL) {

  require(class)
  sitelist <- ppc.fit$sites[[1]]
  
  if(!is.null(nsites)){
    sitelist<-sitelist[1:nsites]
  }
  
  ytr<- data$ytr
  codesy <- levels(ytr)
  
  nc<- dim(split.fit$pr)[3]
  nrows <- trunc(20/nc)
  win.metafile()
  
  par(mfcol=c(nrows,6))
  
  par(mar=c(0,1,0,1))
  par(cex=.6)
  ht <- peak.fit$ind*peak.fit$ht
  
  
  for(i in sitelist){
    
    if (ppc.options$debug) cat(i,fill=T)
    xmax <- max(ht[i,])
    br <- hist(ht[i,],col=3,xlab="",nclass=10,xlim=c(0,xmax*1.15),axes=F,plot=F, main="")$br
    
    par(mar=c(0,1,1,1))
    junk <- hist(ht[i,ytr==codesy[1]],col=3,xlab="",nclass=10,xlim=c(0,xmax*1.15),axes=F,plot=F
                 ,breaks=br)
    hist(ht[i,ytr==codesy[1]],col=5,xlab="",nclass=10,xlim=c(0,xmax*1.15),axes=F,breaks=br, main="")
    box(col=gray)
    wid <- (junk$br[2]-junk$br[1])
    v <- split.fit$cutpoints[i,split.fit$cuthat[i]]+wid
    v0<-v
    if(v==0){v<-v+wid}
    
    h <- max(junk$counts)
    
    abline(v=v,lty=1,col=2)
    
    for(j in (1:ncol(split.fit$prclose))[split.fit$prclose[i,]]){
      points(split.fit$cutpoints[i,j]+wid,h/2,col=2,pch="x")
    }
    
    o<-knn1(data$logmz, centroid.fit$cent[i,1], 1:length(data$logmz))
    lab<- data$mz.labels[o]
    text(xmax*.5,.9*h,labels=lab,cex=.8)
    text(v,3*h/4,labels=as.character(round(v0,2)),col=6)
    
    
    pr <- sum(ht[i,ytr==codesy[1]]>split.fit$cutpoints[i,split.fit$cuthat[i]])/sum(ytr==codesy[1])
    text(xmax,h/2,labels=as.character(round(pr,2)),col=2,cex=.6)
    
    if(nc>2){
      for(ii in 2:(nc-1)){
        par(mar=c(1,1,1,1))
        junk2 <- hist(ht[i,ytr==codesy[ii]],col=4,xlab="",nclass=10,xlim=c(0,xmax*1.15),axes=F,plot=
                      F,breaks=br)
        hist(ht[i,ytr==codesy[ii]],col=ii,xlab="",nclass=10,xlim=c(0,xmax*1.15),axes=F,breaks=br, main="") 
        
        box(col=gray)
        abline(v=v,lty=1,col=2)
        h2 <- max(junk2$counts)
        pr <- sum(ht[i,ytr==codesy[ii]]>split.fit$cutpoints[i,split.fit$cuthat[i]])/sum(ytr==codesy[ii])
        text(xmax,h2/2,labels=as.character(round(pr,2)),col=2,cex=.6)
      }}
    
    par(mar=c(2,1,0,1))
    junk2 <- hist(ht[i,ytr==codesy[nc]],col=4,xlab="",nclass=10,xlim=c(0,xmax*1.15),axes=F,plot=
                  F,breaks=br)
    hist(ht[i,ytr==codesy[nc]],col=nc,xlab="",nclass=10,xlim=c(0,xmax*1.15),axes=F,breaks=br,main="") 
    
    box(col=gray)
    abline(v=v,lty=1,col=2)
    h2 <- max(junk2$counts)
    pr <- sum(ht[i,ytr==codesy[nc]]>split.fit$cutpoints[i,split.fit$cuthat[i]])/sum(ytr==codesy[nc])
    text(xmax,h2/2,labels=as.character(round(pr,2)),col=2,cex=.6)
    
  }
  dev.off()
  return()
}
  
##
## Read a single peaks file and return list of logmz and peak
##
ppc.xl.read.peak.file <- function(filename) {
  ## note: we ignore anything wiht mz <=0 
  pat<-scan(filename,what="", comment.char="#")
  pat<-matrix(pat, ncol=7, byrow=T) [,(1:2)]
  pat<-matrix(as.numeric(as.character(pat)),ncol=2)
  pat <- pat[pat[, 1] > 0, ]
  # rob added log in the next line
  return(cbind(log(pat[,1]),pat[,2]))
}

##
## Read peak data from given files.
##
ppc.xl.read.peak.data <- function(files, batches.exist) {
  ##  if (! batches.exist) {
  n <- length(files)
  peaklist <- vector("list",n)
  mz <- NULL
  for(i  in 1:n){
    w <- ppc.xl.read.peak.file(files[i])
    peaklist[[i]] <- w
    # rob added the exp in the next line
    mz <- c(mz, exp(w[, 1]))
  }
  return(list(peaklist=peaklist, mz=mz))
  ##}
}

##
## Subset and reshape data according to its type
##
ppc.xl.subset.and.reshape <- function(data, user.params) {
  if (data$data.type==ppc.constants$data.types$raw) {
    return(ppc.subset.and.reshape(data, user.params))
  } else {
    return(ppc.subset.and.reshape.peakdata(data, user.params))
  }
}

##
## Depending on the dataset, augment the data with peaks if necessary
##
ppc.xl.augment.with.peaks <- function(data, user.parms) {
  if (data$data.type==ppc.constants$data.types$raw) {
    ## add peaklist to data object
    data$peaklist<- ppc.make.peaklist(data, user.parms)
  }
  return(data)
}


ppc.xl.plotfdr <- function(data, centroid.fit, peak.fit, split.fit, training.fit, user.parms) {
  fdrfit <- ppc.fdr(data, centroid.fit, peak.fit, split.fit, training.fit, user.parms)
  return(list(x=fdrfit$results[,"npeaks"],
              y=fdrfit$results[,"fdr"],
              x.label="No. of significant peaks",
              y.label="False discovery rate"))
}

ppc.xl.remove.suffix <- function(x){
  w <- unlist(strsplit(x, "\\."))
  m <- length(w)
  return(ifelse(m==1, w, paste(w[1:(m-1)], sep="", collapse=".")))
}

