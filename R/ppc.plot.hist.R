ppc.plot.hist <-  function(peak.fit, ppc.fit, centroid.fit, split.fit, data, first.site=1, last.site=25, title.plot=NULL){
# first.site and last.site  indicate how many peak  sites to plot.   default is 
#  1 thru 25

require(class)
sitelist <- ppc.fit$sites[[1]]

if(!is.null(first.site) & is.null(last.site)){
  last.site<- first.site+24
}

if(is.null(first.site) & !is.null(last.site))
             {stop(" Error: first.site missing")}


if(last.site<first.site){stop(" Error: last.site<first.site")}

if(last.site-first.site>24){ last.site<- first.site+24}

if(!is.null(first.site)){
  sitelist<-sitelist[first.site: last.site]
}

ytr<- data$ytr
codesy <- levels(ytr)

nc<- dim(split.fit$pr)[3]
nrows <- trunc(20/nc)
par(mfcol=c(nrows,5))

if(!is.null(title.plot)){ par(oma=c(0,0,5,0))}
par(mar=c(0,1,0,1))
par(cex=.6)
ht <- peak.fit$ind*peak.fit$ht


for(i in sitelist){
 
  if (ppc.options$debug) cat(i,fill=T)
  xmax <- max(ht[i,])

  ## Ensuring that the breaks in each histogram are same for all classes
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
  
  ## Pattern match to figure 
  o<-knn1(data$logmz, centroid.fit$cent[i,1], 1:length(data$logmz))
  lab<- data$mz.labels[o]
  text(xmax*.5,.9*h,labels=lab,cex=.8)
  text(v,3*h/4,labels=as.character(round(v0,2)),col=6)
  
  
  pr <- sum(ht[i,ytr==codesy[1]]>split.fit$cutpoints[i,split.fit$cuthat[i]])/sum(ytr==codesy[1])
  text(xmax,h/2,labels=as.character(round(pr,2)),col=1,cex=.7)
  
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
      text(xmax,h2/2,labels=as.character(round(pr,2)),col=1,cex=.7)
}}

par(mar=c(2,1,0,1))
junk2 <- hist(ht[i,ytr==codesy[nc]],col=4,xlab="",nclass=10,xlim=c(0,xmax*1.15),axes=F,plot=
F,breaks=br)
hist(ht[i,ytr==codesy[nc]],col=nc,xlab="",nclass=10,xlim=c(0,xmax*1.15),axes=F,breaks=br,main="") 
  
box(col=gray)
abline(v=v,lty=1,col=2)
h2 <- max(junk2$counts)
  pr <- sum(ht[i,ytr==codesy[nc]]>split.fit$cutpoints[i,split.fit$cuthat[i]])/sum(ytr==codesy[nc])
text(xmax,h2/2,labels=as.character(round(pr,2)),col=1,cex=.7)

}

if(!is.null(title.plot)){ mtext(title.plot,outer=TRUE,side=3)}

return()
}

