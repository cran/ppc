
ppc.plotfdr <- function(fdrfit){
  plot(fdrfit$results[,"npeaks"],fdrfit$results[,"fdr"],
xlab="Number of peaks called significant",
ylab="False discovery rate",type="b",log="xy")
  axis(3,at=fdrfit$results[,"npeaks"], labels=round(fdrfit$threshold,2))
  
  return()
}
