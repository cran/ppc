ppc.read.peaks.file<- function(filename){
# note: we ignore anything wiht mz <=0 
 pat<-scan(filename,what="", comment.char="#")
        pat<-matrix(pat,ncol=7,byrow=T)[,(1:2)]
        pat<-matrix(as.numeric(as.character(pat)),ncol=2)
        o<-pat[,1]>0
        pat<-pat[o,]
        return( cbind(log(pat[,1]),pat[,2]))
}

