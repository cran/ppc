##
## ppc.cv.error
## Compute error rate, se of error rate and confusion matrices
## from a ppc.cv fit for all values of threshold
##

ppc.cv.error<- function(y, yhat, folds) {
  temp <- ppc.error(y, yhat)
  err <- temp$err
  confusion <- temp$confusion
 
  nfolds <- length(folds)
  n <- length(y)

  codesy <- names(table(y))
  yhat2 <- matrix(codesy[yhat], nrow = nrow(yhat), ncol = ncol(yhat))
  err2 <- matrix(NA, ncol = ncol(yhat), nrow = nfolds)
  temp <- matrix(y, ncol = ncol(yhat), nrow = n)
  ni <- rep(NA, nfolds)

  for(i in 1:nfolds) {
    ii <- folds[[i]]
    ni[i] <- length(folds[[i]])
    err2[i,  ] <- apply(temp[ii,  ] != yhat2[ii,  ], 2, sum) / ni[i]
  }
  se <- sqrt(apply(err2, 2, var) / nfolds)

  return(list(err = err / n,
              confusion = confusion,
              se = se))

}


