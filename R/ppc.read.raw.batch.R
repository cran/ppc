ppc.read.raw.batch<- function(dir, batches, mz=NULL){

comm<- paste("ls ",dir,sep="")

f1<-NULL
for(i in batches){
  temp<- paste(comm, "/",i,sep="")
  temp2<-system(temp,intern=TRUE )
  temp3<-paste(dir,"/",i,"/",temp2,sep="")
  f1<-c(f1,temp3)
}



if(is.null(mz)){
	pat1 <- read.table(f1[1],sep=",")
        mz <- pat1[,1]
}

xtr <- matrix(NA,nrow=length(mz),ncol=length(f1))
for(j in 1:length(f1)){
  if (ppc.options$debug) cat(j,fill=T)
  temp<-  read.table(f1[j],sep=",")

  xtr[,j] <-approx(temp[,1], temp[,2],xout=mz)$y
}

return(list(xtr=xtr,mz=mz,filenames=f1))
}
