
ppc.remove.beforeslash.and.suffix<-
 function(x){
# this function  removes everything before the last "/" at the
#beginning of a filename,  and ".xxx" from the end of a filename
#  of a vector of filenames

for(ii in 1:length(x)){
 for(i in 1:nchar(x[ii])){
     if(substring(x[ii],i,i)=="/") val<-i
 }
 x[ii]<-substring(x[ii],val+1, nchar(x[ii]))
i<-nchar(x[ii])
 while(substring(x[ii],i,i)!="."){
    i<-i-1
     }
 x[ii]<-substring(x[ii],1,i-1)
}
return(x)
}

# Naras's rewrite of Rob's functrion. Doesn't work correctly on linux

##pcc.remove.beforeslash.and.suffix <- function(x) {
##
## this function  removes everything before the last "/" at the
## beginning of a filename,  and ".xxx" from the end of a filename
## of a vector of filenames
##
##base.names <- basename(x)
##return(sapply(base.names, ppc.remove.suffix))
##}


