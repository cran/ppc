##
## This function min.na is not needed as the function min offers a direct
## facility for not considering NAs.
##
##min.na<-function(x) {
##  min(x[!is.na(x)])
##}

which.is.min <- function(x) {
  y <- seq(length(x))[x == min(x, na.rm=TRUE)]
  y <- y[!is.na(y)]
  if(length(y) > 1)
    y <- sample(y, 1)
  y
}

medoid <- function(x) {
  n <- length(x)
  if(n==1){med <- x}
  if(n>1){
    d <- dist(x)
    dd <- matrix(0,nrow=n,ncol=n)
    dd[row(dd)>col(dd)] <- d
    dd <- dd+t(dd)
    
    m <- apply(dd,2,sum)
    med <- x[(1:n)[m==min(m)]]
    if(length(med)>1){
      med <- mean(med)
    }
  }
  return(med)
}

permute.rows <-function(x) {
  dd <- dim(x)
  n <- dd[1]
  p <- dd[2]
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}


balanced.folds <- function(y, nfolds = min(min(table(y)), 10)) {
  totals <- table(y)
  fmax <- max(totals)
  nfolds <- min(nfolds, fmax)     
  ## makes no sense to have more folds than the max class size
  folds <- as.list(seq(nfolds))
  yids <- split(seq(y), y)        
  ## nice we to get the ids in a list, split by class
  ##Make a big matrix, with enough rows to get in all the folds per class
  bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
  for(i in seq(totals)) {
    bigmat[seq(totals[i]), i] <- sample(yids[[i]])
  }
  smallmat <- matrix(bigmat, nrow = nfolds)       # reshape the matrix
  ## Now do a clever sort to mix up the NAs
  smallmat <- permute.rows(t(smallmat))
  ## Now a clever unlisting
  ## the "clever" unlist doesn't work when there are no NAs
  ##       apply(smallmat, 2, function(x)
  ##        x[!is.na(x)])
  res <-vector("list", nfolds)
  for(j in 1:nfolds) {
    jj <- !is.na(smallmat[, j])
    res[[j]] <- smallmat[jj, j]
  }
  return(res)
}

ppc.peaks <- function(x,span){
  
  
# note-  changed this so that span is now a percentage rather than a
#   number of points
 
   ispan<-trunc(length(x)*span)

  if(ispan%%2==0){ispan <- ispan+1}
  
  n <- length(x) 
  
  junk <- .Fortran("peaks",
                   x,
                   as.integer(ispan),
                   as.integer(n),
                   ans=integer(n),
                   PACKAGE="ppc")
  return(junk$ans==1)
}

hclust.1d <- function(x, debug=FALSE) {
  ##(fast) complete linkage hierarhical clustering, in one dimension
  ## R. Tibshirani Nov 2003
  # to make it work under Linux,
  
  n<-length(x)
  storage.mode(x)<- "double"
  storage.mode(n)<- "integer"
  ##
  ##

  idebug <- ifelse(debug, 1, 0)
  
  junk<- .Fortran("hclust1d",
                  as.integer(idebug),
                  x,
                  n,
                  merge=integer( (n-1)*2),
                  height=double(n-1),
                  order=integer(n),
                  scrat=double(n*3),
                  scrat2=double(n*3),
                  PACKAGE="ppc")

  
  merge<-matrix(junk$merge,ncol=2,byrow=F)
  order<-junk$order
  height<-junk$height
  return(list(merge=merge,height=height,order=order))
}





