ppc.read.peaks.batch<- function(dir, batches){

comm<- paste("ls ",dir,sep="")

f1<-NULL
for(i in batches){
  temp<- paste(comm, "/",i,sep="")
  temp2<-system(temp,intern=TRUE )
  temp3<-paste(dir,"/",i,"/",temp2,sep="")
  f1<-c(f1,temp3)
}


n<-length(f1)
peaklist<- vector("list",n)

logmz<-NULL
 for(i  in 1:length(f1)){
        peaklist[[i]]<- ppc.read.peaks.file(f1[i])
        logmz<-c(logmz,peaklist[[i]][,1])
}

return(list(peaklist=peaklist,filenames=f1, logmz=logmz))

}
