ppc.find.splits<- function(centroid.fit, peak.fit, data, user.parms)
{
  ## find best split points for training data
  ## takes centroids.fit- result of call to make.centroids.list
  ## and peaks.fit- result of call to predict.peaks
  ##
  ##  ytr is  the vector of  class labels 
  ## nsplits is number of equally spaced  split points to try (plus the value 0)
  ##  fix.at.one=-TRUE means zero is the only split point tried
  ##    (ie no peak vs peak)
  ##
  ## returns prhat - estimated optimal split proportions for each site, 
  ##  pr- proportions in each  split category for each class, 
  ## prclose- logical matrix indicating split values that coem within 10% of best
  ## cuthat,cutpoints- optimal cut ## and matrix of all cutpoints considered
  ##
  
  ytr <- data$ytr
  
  n.class<- table(ytr)
  
  nsplits <- user.parms$nsplits
  fix.at.one <- user.parms$fix.at.one 
  
  which.is.max.na <- function(x)
    {
      xx <- x[!is.na(x)]
      y <- seq(length(xx))[xx == max(xx)]
      if(length(y) > 1) {
        y <- sample(y, 1)
      }
      o <- (1:length(x))[!is.na(x)]
      return(o[y])
    }
  ht <- peak.fit$ht * peak.fit$ind
  cent <- centroid.fit$cent
  alpha <- (1:nsplits)/(nsplits + 1)
  p <- nrow(ht)
  Y <- model.matrix( ~ factor(ytr) - 1, data = list(ytr = ytr))
  K <- length(table(ytr))
  pr <- array(0, c(p, nsplits + 1, K))
  qu <- matrix(NA, nrow = p, ncol = nsplits + 1)
  
  ##at each site, compute proportion  of peak heights exceeding the split  quantiles 
  
  for(i in 1:p) {
    if (ppc.options$debug) cat(i,fill=T)
    o <- ht[i,  ] > 0
    temp <- sort(ht[i, o])
    temp2 <- pmax(1, trunc(length(temp) * alpha))
    qu[i,  ] <- c(0, temp[temp2])
    for(j in 1:(nsplits + 1)) {
      for(k in 1:K) {
        pr[i, j, k] <- sum(ht[i, Y[, k] == 1] > qu[i, j])/sum(Y[, k] == 1)
      }
    }
  }
  prmean <- apply(pr, c(1, 2), mean)
  
  ## find best split indices
  
  if(!fix.at.one) {
    prtemp <- apply(abs(pr - array(prmean, c(p, nsplits + 1, K))),
                    c(1, 2), sum)
    cuthat <- apply(prtemp, 1, which.is.max.na)
  }
  else {
    cuthat <- rep(1, p)
  }
  prhat <- matrix(NA, nrow = p, ncol = K)
  
  
  ##compute optimized proportions and  splits within 10% of the best
  
  if(fix.at.one){prclose <- NULL}
  
  
  prclose <- matrix(F, nrow = p, ncol = nsplits + 1)
  for(i in 1:p) {
    prhat[i,  ] <- pr[i, cuthat[i],  ]
    if(!fix.at.one){
      for(j in 1:(nsplits + 1)) {
        prclose[i, j] <- prtemp[i, j] > 0.9 * prtemp[i, cuthat[
                                                               i]]
      }}
  }
  
  return(list(prhat = prhat, pr = pr, n.class=n.class, cutpoints = qu, cuthat = cuthat, prclose
              = prclose,nsplits=nsplits, fix.at.one=fix.at.one))
}


