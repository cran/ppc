ppc.predict <- function(centroid.fit, split.fit, logmz,  peaklist.te, n.threshold = 30, threshold=NULL, metric = 
                        c("binomial","euclidean", "absolute"), summ=c("mean","median"))
{
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

SMALL<- 10e-6

  metric <- match.arg(metric)
  print(metric)
summ <- match.arg(summ)
  print(summ)
  if(!is.null(threshold)){n.thresholds<-length(threshold)}
  
  n <- length(peaklist.te)
  n.class <- split.fit$n.class
  K<- ncol(split.fit$prhat)
  p<-length(split.fit$cuthat)
  yhat <- matrix(NA, nrow = p, ncol = n)
  yhatt <- matrix(NA, ncol = n.threshold, nrow = n)
  ind0 <- yhat
  ht <- yhat
  dis <- array(NA,c(n, K, n.threshold))
  
  numsites <- rep(NA, n.threshold)
  sites <- vector("list", n.threshold)
  for(j in 1:n) {
    if (ppc.options$debug) cat(j)
    aaa <- ppc.predict.peaks1(centroid.fit, logmz,  peaklist= 
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
  
  if(is.null(threshold)){threshold <- seq(0, max(abs(delta)), length = n.threshold)}
  
  soft.thresh <- function(x, tt)
    {
      sign(x) * (abs(x) - tt) * (abs(x) > tt)
    }
  
  for(j in 1:n.threshold) {
    if (ppc.options$debug) cat(j)
    delta2 <- soft.thresh(delta, threshold[j])
    pr <- prmean + delta2
    sumabs<-  apply(abs(delta2),1,sum)
    pos <-  sumabs> SMALL
    numsites[j] <- sum(pos)
    sites[[j]] <- (1:nrow(pr))[pos]
    
    sites[[j]] <- sites[[j]][order( - abs(sumabs[pos]))]
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
    
    yhatt[, j] <- apply(dd,1,which.is.min)
    dis[,,j] <- dd
    
  }
  prob <- exp(-dis)
  prob <- prob/array(apply(prob,c(1,3),sum),dim(prob))
  
  return(list(yhat = yhatt, threshold = threshold, numsites = numsites, sites
              = sites, ind = yhat, ind0=ind0, prob=prob, ht=ht))
  
}
