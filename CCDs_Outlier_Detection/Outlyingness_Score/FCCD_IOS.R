source("/mmfs1/home/rzs0112/code_working_folder/Outlyingness_Score/Outlyingness_Score.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/RK_CCD_New.R")

# calculate the Inbound Outlyingness Scores of RKCCD
# niter is the number of simulated data set for CSR test if simul is not provided
# d: dimensionality
# niter: the number of data sets to simulate if simul is null
FCCD_IOS = function(datax, simul=NULL, d, quant=0.99){
  # The digraph based on RKCCD
  Mgraph = rccd.clustering.mutual.connected_correct(datax, low.num=2, r.seq=10, 
                                                    dom.method="greedy2", quan=quant, simul=simul, niter=1000, scores=T)
  # find the optimal number of clusters by maximizing silhouette index
  result_cls = rccd.silhouette_mutual(Mgraph,datax,ind=NULL, lenClimit=Inf,
                                      k=max(d,5), min.cls = max(0.05,cont/2))
  # The radius of each point
  R = Mgraph$R[order(Mgraph$D)]
  vd = Vic_Den(datax,R,d) # the vicinity density
  # Number of clusters
  n_cls = result_cls$sil.ind
  member = unique(result_cls$label)
  # calculate the scores for each cluster and standardization
  scores = lapply(member,function(x){
    index_cls = which(result_cls$label==x)
    cluster = datax[index_cls,]
    R_cls = R[index_cls]
    score_cls = IOS(cluster,R_cls,d)
    score_cls = std_MADN(score_cls)
    return(score_cls)
  })
  scores_whole = rep(0,length(datax[,1]))
  for(i in 1:n_cls){
    scores_whole[which(result_cls$label==member[i])] = scores[[i]]
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