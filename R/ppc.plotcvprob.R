ppc.plotcvprob <- function(fit, data, threshold) {
  par(pch = 1)
  ii <- (1:length(fit$threshold))[fit$threshold > threshold]
  ii <- ii[1]
  ss <- data$samplelabels
  pp <- fit$prob[,  , ii]

  y <- fit$y

  o <- order(y)
  y <- y[o]
  if(!is.null(ss)) {
    ss <- ss[o]
  }
  ppp <- pp[o,  ]
  n <- nrow(ppp)
  nc <- length(unique(y))
  par(cex = 1)
  plot(1:n, ppp[, 2], type = "n", xlab = "sample", ylab = 
       "cross-validated probabilities",  axes = FALSE)
  axis(1)
  labs <- round(seq(min(ppp),max(ppp),length=5),2)
  axis(2, labels = as.character(labs),at=labs)
  for(j in 1:nc) {
    points(1:n, ppp[, j], col = j + 1)
  }
  for(j in 1:(nc - 1)) {
    abline(v = cumsum(table(y))[j] + 0.5, lty = 2)
  }
  h <- c(0, table(y))
  for(j in 2:(nc + 1)) {
    text(sum(h[1:(j - 1)]) + 0.5 * h[j], 1.02, label = levels(y)[j - 
                                                 1], col = j)
  }
  abline(h = 1)
  if(!is.null(ss)) {
    text(1:length(ss), 1.1, labels = ss, srt = 90, cex = 0.7)
  }
  ##if(!is.null(ss)){axis(3,labels=ss,at=1:length(ss),srt=90)}
}

