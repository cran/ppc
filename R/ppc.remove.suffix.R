

ppc.remove.suffix<-
 function(x){
# this function just removes a .xxx suffix from the end of a filename
#  of a vector of filenames

for(ii in 1:length(x)){
 i<-nchar(x[ii])
while(substring(x[ii],i,i)!= "."){i<-i-1}
 x[ii]<-substring(x[ii],1,i-1)
}
return(x)
}



