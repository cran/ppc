ppc.make.centroid.list<-
function (data, user.parms) 
{
    peaklist <- data$peaklist
    ytr <- data$ytr
    logmz <- data$logmz
    recluster <- user.parms$recluster
    peak.gap <- user.parms$peak.gap
    n <- length(ytr)
    m <- NULL
    for (i in 1:n) {
        m <- c(m, peaklist[[i]][, 1])
    }
    mm <- sort(unique(m))
    clust.tree <- hclust.1d(mm, debug = as.integer(ppc.options$debug))
    if (ppc.options$debug) 
        cat("bef ppc.make.centroids", fill = TRUE)
    cent <- ppc.make.centroids(clust.tree, mm, user.parms)
    if (ppc.options$debug) 
        cat("aft ppc.make.centroids", fill = TRUE)
    return(list(cent = cent, peaklist = peaklist, clust.tree = clust.tree, 
        all.peaks = mm, peak.gap = user.parms$peak.gap, recluster = user.parms$recluster))
}
