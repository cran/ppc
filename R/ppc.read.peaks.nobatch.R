ppc.read.peaks.nobatch<- function(dir){
comm<- paste("ls ",dir,sep="")
f1<-system(comm,intern=T)

n<-length(f1)
peaklist<- vector("list",n)

logmz<-NULL
 for(i  in 1:length(f1)){
        peaklist[[i]]<- ppc.read.peaks.file(paste(dir, "/",f1[i],sep=""))
  logmz<-c(logmz,peaklist[[i]][,1])
}

return(list(peaklist=peaklist,filenames=f1, logmz=logmz))

}
