 ppc.find.peaks  <- function(mz,x,user.parms){
##
## peak finding algorithm for SELDI spectra, from 
##  description of Yasui  et al Biostatistics 2003
##
## x is a single spectrum
##  span is window size for supersmoother estimate of background
##         minht is min height for a peak; stn is  min signal to noise
##         ratio for a peak ( i.e. must be > stn* background)
  
  span<-user.parms$span
  minht<-user.parms$minht
  
  stn<-user.parms$stn
  smoothing.span<-user.parms$smoothing.span
  
  a<-ppc.peaks(x,span=span)
  b<-supsmu(mz,x,span=smoothing.span)
  peaks2<-(1:length(mz))[a & (x> stn*b$y) & (x>minht)]
  return (cbind(mz[peaks2],x[peaks2]))
}
