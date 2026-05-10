source("/mmfs1/home/rzs0112/code_working_folder/Outlyingness_Score/Outlyingness_Score.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/UN_CCD.R")

# calculate the Outbound Outlyingness Scores of NNCCD
# niter is the number of simulated data set for CSR test if simul is not provided
# method, accend: select the radius in accend order, decent: select the radius in decent order
# d: dimensionality
NNCCD_OOS = function(datax, simul=NULL, method="ascend",d){
  # The digraph based on NNCCD
  graph <- nnccd_clustering_quantile(datax, low.num=3, quantile="lower", 
                                     method=method, dom.method="greedy2", 
                                     simul=simul, niter=1000, scores=T)
  # The radius of each point
  R = graph$R[order(graph$D)]
  
  # Calculating the OSS of the entire data set and standardization
  scores = OOS(datax,R,d)
  scores = std_MADN(scores)
  vd = Vic_Den(datax,R,d) # the vicinity density
  
  # break ties
  score_vd = cbind(scores,vd)
  order_index = order(scores)
  score_vd = score_vd[order_index,]
  frequency = table(scores)
  repeated_values = as.numeric(names(frequency[frequency > 1]))
  for(x in repeated_values){
    if(x==0){
      index = which(score_vd[,1]==x)
    } else {
      index = which(abs(score_vd[,1]-x)/x<0.0001)
    }
    if(length(index)<2) break
    index_max = max(index)
    index_min = min(index)
    if(index_min==1){
      s = score_vd[,1][1]
      e = score_vd[,1][index_max+1]
      diff = e - s
    } else if(index_max==length(datax[,1])){
      s = score_vd[,1][index_min-1]
      e = score_vd[,1][index_max]
      diff = e - s
    } else {
      s = score_vd[,1][index_min-1]
      e = score_vd[,1][index_max+1]
      diff = e - s
    }
    if(is.na(diff)) break
    den_sum = sum(score_vd[,2][index])
    score_vd[,1][index] = sapply(index, function(x){
      new_s = s + diff*score_vd[,2][x]/den_sum
      return(new_s)
    })
    scores = score_vd[order(order_index),1]
  }
  return(scores)
}

# calculate the Outbound Outlyingness Scores of NNCCD
# niter is the number of simulated data set for CSR test if simul is not provided
# method, ascend: select the radius in accend order, descend: select the radius in decent order
# d: dimensionality
NNCCD_IOS = function(datax, simul=NULL, method="ascend", d, min.cls=0){
  # The digraph based on NNCCD
  graph <- nnccd_clustering_quantile(datax, low.num=3, quantile="lower", 
                                     method=method, dom.method="greedy2", 
                                     simul=simul, niter=1000, scores=T)
  # The radius of each point
  R = graph$R[order(graph$D)]
  # The clustering result and labeling
  NN_cluster = nnccd.silhouette(graph,datax,cls=NULL, min.cls=min.cls, ind=NULL, lenDlimit=Inf)
  vd = Vic_Den(datax,R,d) # the vicinity density
  # Number of clusters
  n_cls = NN_cluster$si.ind
  member = unique(NN_cluster$label)
  # calculate the scores for each cluster and standardization
  scores = lapply(member,function(x){
    index_cls = which(NN_cluster$label==x)
    cluster = datax[index_cls,]
    R_cls = R[index_cls]
    score_cls = IOS(cluster,R_cls,d)
    score_cls = std_MADN(score_cls)
    return(score_cls)
  })
  scores_whole = rep(0,length(datax[,1]))
  for(i in 1:n_cls){
    scores_whole[which(NN_cluster$label==member[i])] = scores[[i]]
  }
  
  # break ties
  score_vd = cbind(scores_whole,vd)
  order_index = order(scores_whole)
  score_vd = score_vd[order_index,]
  frequency = table(scores_whole)
  repeated_values = as.numeric(names(frequency[frequency > 1]))
  for(x in repeated_values){
    if(x==0){
      index = which(score_vd[,1]==x)
    } else {
      index = which(abs(score_vd[,1]-x)/x<0.0001)
    }
    if(length(index)<2) break
    index_max = max(index)
    index_min = min(index)
    if(index_min==1){
      s = score_vd[,1][1]
      e = score_vd[,1][index_max+1]
      diff = e - s
    } else if(index_max==length(datax[,1])){
      s = score_vd[,1][index_min-1]
      e = score_vd[,1][index_max]
      diff = e - s
    } else {
      s = score_vd[,1][index_min-1]
      e = score_vd[,1][index_max+1]
      diff = e - s
    }
    if(is.na(diff)) break
    den_sum = sum(score_vd[,2][index])
    score_vd[,1][index] = sapply(index, function(x){
      new_s = s + diff*score_vd[,2][x]/den_sum
      return(new_s)
    })
    scores_whole = score_vd[order(order_index),1]
  }
  return(scores_whole)
}