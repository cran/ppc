ppc.plotcv <- function(fit) {
  par(mar = c(5, 5, 5, 1))
  par(mfrow = c(2, 1))
  n <- nrow(fit$yhat)
  y <- fit$y
  nc <- length(table(y))
  nfolds <- length(fit$folds)

  plot(fit$threshold, fit$err, ylim = c(-0.1, 0.8), xlab = 
       "Value of threshold  ", ylab = "Misclassification Error", type
       = "n", yaxt = "n")
  axis(3, at = fit$threshold, lab = paste(fit$numsites), srt = 90, adj = 0)
  mtext("Number of sites", 3, 4, cex = 1.2)
  axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8))
  lines(fit$threshold, fit$err, col = 2)
  o <- fit$err == min(fit$err)
  points(fit$threshold[o], fit$err[o], pch = "x")
  error.bars(fit$threshold, fit$err - fit$se, fit$err + fit$se)
  err2 <- matrix(NA, nrow = length(unique(y)), ncol = length(fit$threshold
                                                 ))
  for(i in 1:(length(fit$threshold) - 1)) {
    s <- fit$confusion[[i]]
    diag(s) <- 0
    err2[, i] <- apply(s, 1, sum)/table(y)
  }
  plot(fit$threshold, err2[1,  ], ylim = c(-0.1, 1.1), xlab = 
       "Value of threshold ", ylab = "Misclassification Error", type
       = "n", yaxt = "n")
  axis(3, at = fit$threshold, lab = paste(fit$numsites), srt = 90, adj = 0)   
 # mtext("Number of sites", 3, 4,cex=1.2)
  axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8))
  for(i in 1:nrow(err2)) {
    lines(fit$threshold, err2[i,  ], col = i + 1)
  }
  legend(0, 0.9, dimnames(table(y))[[1]], col = (2:(nc + 1)), lty = 1)
  par(mfrow = c(1, 1))
}

error.bars <-function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}

