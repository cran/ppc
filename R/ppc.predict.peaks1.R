ppc.predict.peaks1 <-  function(centroid.fit, logmz,  peaklist.new)
{
  require(class)
  ## looks in peaklist.new for peaks that are near the centroid peaks
  ##   defined in centroid.fit
  ## returns logical indicator vector "ind" and peak heights "ht"
  
# changed to work with new test data with new mz values
#  pee <- match(peaklist.new[, 1], logmz)

pee<-knn1(matrix(logmz,ncol=1),peaklist.new[, 1,drop=FALSE], 1:length(logmz))
  
  xnew <- rep(NA, length(logmz))
  xnew[pee] <- peaklist.new[, 2]
  
  m3 <- logmz[pee]
  cent <- centroid.fit$cent
 
  bb <- knn1(matrix(m3, ncol = 1), cent[, 1, drop = F], 1:length(m3))
  
  ind <- abs(m3[bb] - cent[, 1]) <=centroid.fit$peak.gap/2 
  ht <- xnew[pee][bb]

#NOTE: ht is returned as non-zero, even if ind=0 (ie peak does not match
#a training peak)

  return(list( ind =ind, ht = ht))
}
