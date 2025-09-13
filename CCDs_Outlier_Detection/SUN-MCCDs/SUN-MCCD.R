source("/mmfs1/home/rzs0112/code_working_folder/ccds/UN_CCD.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/mKNN_CCD_functions.R")


# datax: the data set for analysis
# Simul:  Simulated boundaries for K-test
# cont: (expected) percentage of clusters
# low.num: the least number of points in each covering ball
SUNMCCD_outlier = function(datax, simul=NULL, min.cls=0, method="ascend", low.num=3){
  
  graph <- nnccd_clustering_quantile(datax, low.num=low.num, quantile="lower", 
                                     method=method, dom.method="greedy2", 
                                     simul=simul, niter=100)
  NN_result = nnccd.silhouette(graph, datax, cls=NULL, min.cls=min.cls, ind=NULL, lenDlimit=Inf)
  NN_result = c(graph,NN_result)
  
  # find the covering balls that are mutually catched with the dominating covering balls
  r = NN_result$R[order(NN_result$D)]
  dom_index = NN_result$Int.D[1:NN_result$si.ind] # the index of the dominating covering balls
  ddatax = as.matrix(dist(datax))
  n = dim(datax)[1]
  
  M <- matrix(as.integer(ddatax <= r), length(r))
  dom_index_inter=lapply(dom_index,function(i){
    temp_list = sapply(1:n,function(j){
      if(M[i,j] & M[j,i])
        return(j)
    })
    return(unlist(temp_list))
  })
  
  # a list of clusters of the dataset
  nnccd_clusters = lapply(1:NN_result$si.ind, function(x){return(datax[NN_result$label==x,])})
  # the cover info of each cluster
  catch_index = lapply(1:NN_result$si.ind, function(x){
    index_temp = dom_index_inter[[x]]
    temp_list = sapply(1:n,function(j){
      if(any(M[index_temp,j]))
        return(j)
    })
    return(unlist(temp_list))
  })
  
  catch_points = lapply(1:NN_result$si.ind, function(x){datax[catch_index[[x]],]})
  
  #find the max density delta such that the KS-CCD for each Dominated graph is connected
  delta_max = sapply(catch_points, connected.ksccd.m)
  
  # find the connected components for each cluster with the delta obtained above.
  cluster_component = lapply(1:NN_result$si.ind, function(x){
    return(ksccd.connected(nnccd_clusters[[x]],delta_max[x],sequential=FALSE,alpha=0.05)$member)})
  return(c(clusters=list(nnccd_clusters),label=list(cluster_component)))
}