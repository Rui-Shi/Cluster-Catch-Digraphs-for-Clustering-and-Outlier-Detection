source("/mmfs1/home/rzs0112/code_working_folder/ccds/UN_CCD.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/mKNN_CCD_functions.R")

# main function for outlier function
# datax: the data set
# simul: the simulated windows or envelope 
UNMCCD_outlier = function(datax, simul, method="ascend", niter=1000){
  graph <- nnccd_clustering_quantile(datax, low.num=3, quantile="lower", 
                                     method=method, dom.method="greedy2", 
                                     simul=simul, niter=niter)
  NN.result = nnccd.silhouette(graph,datax,cls=NULL,ind=NULL, lenDlimit=Inf)
  NN.result = c(graph,NN.result)
  
  # a list of clusters of the dataset
  nnccd.clusters = lapply(1:NN.result$si.ind, function(x){return(datax[NN.result$label==x,])})
  
  #the cover info of each Dominated graph
  catch.num = NN.result$catch
  ddatax = as.matrix(dist(datax))
  catch.info = lapply(1:NN.result$si.ind, function(x){return(datax[ddatax[NN.result$Int.D[x],]<=NN.result$Int.R[x],])})
  
  #find the max density delta such that the KS-CCD for each Dominated graph is connected
  delta.max = sapply(catch.info, connected.ksccd.m)
  
  # find the connected components for each cluster with the delta obtained above.
  cluster.component = lapply(1:NN.result$si.ind, function(x){
    return(ksccd.connected(nnccd.clusters[[x]],delta.max[x],sequential=FALSE,alpha=0.05)$member)})
  return(c(clusters=list(nnccd.clusters),label=list(cluster.component)))
}