source("/mmfs1/home/rzs0112/code_working_folder/ccds/RK_CCD_New.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/mKNN_CCD_functions.R")


# datax: the data set for analysis
# Simul:  Simulated boundaries for K-test
# cont: (expected) percentage of clusters
# low.num: the least number of points in each covering ball
SUMCCD_outlier = function(datax, simul=NULL, min.cls=0, low.num=2, quant=0.99){
  RK_result = RKCCD_correct_quant(datax, low.num=low.num, r.seq=10, dom.method="greedy2", 
                                  quan=quant, simul=simul, cls=NULL,ind=NULL, niter=100, min.cls=min.cls)
  
  # find the covering balls that are mutually catched with the dominating covering balls
  r = RK_result$R[order(RK_result$D)]
  dom_index = RK_result$Int.D[1:RK_result$si.ind] # the index of the dominating covering balls
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
  rkccd_clusters = lapply(1:RK_result$si.ind, function(x){return(datax[RK_result$label==x,])})
  # the cover info of each cluster
  catch_index = lapply(1:RK_result$si.ind, function(x){
    index_temp = dom_index_inter[[x]]
    temp_list = sapply(1:n,function(j){
      if(any(M[index_temp,j]))
        return(j)
    })
    return(unlist(temp_list))
  })
  
  catch_points = lapply(1:RK_result$si.ind, function(x){datax[catch_index[[x]],]})
  
  #find the max density delta such that the KS-CCD for each Dominated graph is connected
  delta_max = sapply(catch_points, connected.ksccd.m)
  
  # find the connected components for each cluster with the delta obtained above.
  cluster_component = lapply(1:RK_result$si.ind, function(x){
    return(ksccd.connected(rkccd_clusters[[x]],delta_max[x],sequential=FALSE,alpha=0.05)$member)})
  return(c(clusters=list(rkccd_clusters),label=list(cluster_component)))
}