ppc.make.centroids <-function(clust.tree, x, user.parms)
{
  peak.gap <-  user.parms$peak.gap
  nclust <- user.parms$nclust
  recluster <- user.parms$recluster
  
  if(!is.null(peak.gap) & !is.null(nclust))
    stop("Error: only one of peak.gap  and nlcust should be specified"
         )
  if(!is.null(nclust)) {
    aa <- cutree(clust.tree, k = nclust)
  }
  if(!is.null(peak.gap)) {
    aa <- cutree(clust.tree, h = peak.gap)
    n <- nrow(clust.tree$merge) + 1
    uaa <- unique(aa)
    cen <- matrix(NA, ncol = 4, nrow = length(uaa))
    i <- 0
    for(ii in uaa) {
      b <- x[(1:n)[aa == ii]]
      i <- i + 1
                                        #     cen[i, 1] <- (min(b) + max(b))/2
                                        #      cen[i, 2] <- (max(b) - min(b))/2
      if (ppc.options$debug) cat(i, fill = T)
      cen[i, 1] <- medoid(b)
                                        # i changed this on oct 27
                                        # cen[i, 2] <- max(cen[i,1] - min(b), max(b) - cen[i,1])
      cen[i, 2] <- max(b) - min(b)
      cen[i, 3] <- min(b)
      cen[i,4] <-  max(b)
    }
  }
  
  if(recluster & !is.null(peak.gap)){
    mind<- min(diff(cen[,1]))
    while(mind < peak.gap){
      if (ppc.options$debug) cat("reclustering",fill=T)
      xx <- cen
      clust.tree<- hclust.1d(xx[,1], debug=ppc.options$debug)
      aa <- cutree(clust.tree, h = peak.gap)
      n <- nrow(clust.tree$merge) + 1
      uaa <- unique(aa)
      cen <- matrix(NA, ncol = 4, nrow = length(uaa))
      i <- 0
      for(ii in uaa) {
        b <- xx[,1][(1:n)[aa == ii]]
        bmin <- xx[,3][(1:n)[aa == ii]]
        bmax <- xx[,4][(1:n)[aa == ii]]
        i <- i+1
        
        cen[i, 1] <- mean(b)
        cen[i, 3] <- min(bmin)
        cen[i,4] <-  max(bmax)
        cen[i,2] <- cen[i,4]-cen[i,3]
      }
      mind<- min(diff(cen[,1]))
    }
  }
  
  dimnames(cen) <- list(NULL, c("pos", "diameter", "min", "max"))
  o <- order(cen[, 1])
  return(cen[o,  ])
}

