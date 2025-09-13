source("/media/rui/exNVME/code_working_folder/ccds/RK_CCD_New.R")
source("/media/rui/exNVME/code_working_folder/ccds/mKNN_CCD_functions.R")

# datax: the data set for analysis
# Simul:  Simulated boundaries for K-test
# cont: (expected) percentage of clusters
# low.num: the least number of points in each covering ball
# k: the number to neighbors to compute LDOF
MFCCD_outlier = function(datax, simul=NULL, min.cls=0, low.num=2, quant=0.99, k=NULL){
  # datax=data.list[[5]]
  # min.cls=0.04; low.num=2; quant=0.99; k=NULL
  n = dim(datax)[1]
  d = dim(datax)[2]
  if(is.null(k)){k=max(d,5)}
  RK_result = RKCCD_correct_quant(datax, low.num=low.num, r.seq=10, dom.method="greedy2", 
                                  quan=quant, simul=simul, cls=NULL,ind=NULL, niter=100, min.cls=min.cls)
  
  # find the covering balls that are mutually catched with the dominating covering balls
  r = RK_result$R[order(RK_result$D)]
  dom_index = RK_result$Int.D # the index of the dominating covering balls
  ddatax = as.matrix(dist(datax)) # dist matrix
  dist.M = ddatax
  
  M <- matrix(as.integer(ddatax <= r), length(r))
  # find the index of the points that are mutual catched with dominating covering balls
  dom_index_inter=lapply(dom_index,function(i){
    temp_list = sapply(1:n,function(j){
      if(M[i,j] & M[j,i])
        return(j)
    })
    return(unlist(temp_list))
  })
  
  # a list of clusters of the dataset
  # the cover info of each union of covering balls
  catch_index = lapply(1:length(dom_index_inter), function(x){
    index_temp = dom_index_inter[[x]]
    temp_list = sapply(1:n,function(j){
      if(any(M[index_temp,j]))
        return(j)
    })
    return(unlist(temp_list))
  })
  
  # find the points that are covered more than once
  dup_index = unique(unlist(catch_index)[duplicated(unlist(catch_index))])
  
  # member is the labels for the data set
  member = rep(0,n)
  for(i in 1:length(catch_index)) member[catch_index[[i]]]=i
  member[dup_index]=0
  
  catch_points = lapply(1:length(dom_index_inter), function(x){return(datax[catch_index[[x]],])})
  counts = table(member[which(member>0)])
  D.member = as.numeric(names(counts)[order(counts,decreasing = T)]) # the labels in decenting order
  graph = list(member=member, D.member=D.member, dist.M=ddatax)
  
  # assign unlabeled points with LDOF and update memberships, resulting graph_full
  graph_full = graph
  graph_full$member=rccd.clustering.nonvalid.knn1(graph, datax, k)
  counts = table(graph_full$member[which(graph_full$member>0)])
  graph_full$D.member = as.numeric(names(counts)[order(counts,decreasing = T)])
  # Update the labels of duplicate point
  graph$member[dup_index]=graph_full$member[dup_index]
  counts = table(graph$member[which(graph$member>0)])
  graph$D.member=as.numeric(names(counts)[order(counts,decreasing = T)])
  
  # find optimal clustering by max average Silhouette index
  result_cls = rccd.silhouette_mutual1(graph=graph_full, datax=datax, k=k, min.cls = min.cls)
  
  
  #find the max density delta such that the KS-MCG for covered cluster is connected
  catch_index = lapply(1:result_cls$sil.ind, function(x){
    label_temp = graph_full$D.member[x]
    catch_index_temp = which(graph$member==label_temp)
    return(catch_index_temp)
  })
  catch_points = lapply(1:length(catch_index), function(x){return(datax[catch_index[[x]],])})
  
  delta_max = sapply(catch_points, connected.ksccd.m)
  
  # find the connected components for each cluster with the delta obtained above.
  cluster_full_index = lapply(1:result_cls$sil.ind, function(x){
    label_temp = graph_full$D.member[x]
    cluster_index_temp = which(result_cls$label==label_temp)
    return(cluster_index_temp)
  })
  cluster_full = lapply(1:length(cluster_full_index), function(x){return(datax[cluster_full_index[[x]],])})
  
  cluster_component = lapply(1:result_cls$sil.ind, function(x){
    return(ksccd.connected(cluster_full[[x]],delta_max[x],sequential=FALSE,alpha=0.05)$member)})
  return(c(clusters=list(cluster_full),label=list(cluster_component)))
}