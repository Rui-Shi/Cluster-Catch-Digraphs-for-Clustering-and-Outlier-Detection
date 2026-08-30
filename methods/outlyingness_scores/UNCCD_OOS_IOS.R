source(here::here("methods/outlyingness_scores/Outlyingness_Score.R"))
source(here::here("R/ccds/UN_CCD.R"))

# calculate the Outbound Outlyingness Scores of NNCCD
# niter is the number of simulated data set for CSR test if simul is not provided
# method, accend: select the radius in accend order, decent: select the radius in decent order
# d: dimensionality
NNCCD_OOS = function(datax, simul=NULL, method="ascend",d, tiebreak=getOption("ccd.tiebreak","eq10"), parts=FALSE){
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
  raw_scores = scores
  scores = break_ties(scores, vd, mode = tiebreak)
  if(parts){
    return(list(score=scores, raw=raw_scores, vd=vd, label=NULL))
  }
  return(scores)
}

# calculate the Outbound Outlyingness Scores of NNCCD
# niter is the number of simulated data set for CSR test if simul is not provided
# method, ascend: select the radius in accend order, descend: select the radius in decent order
# d: dimensionality
NNCCD_IOS = function(datax, simul=NULL, method="ascend", d, min.cls=0, tiebreak=getOption("ccd.tiebreak","eq10"), parts=FALSE){
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
  raw_scores = scores_whole
  scores_whole = break_ties(scores_whole, vd, mode = tiebreak)
  if(parts){
    return(list(score=scores_whole, raw=raw_scores, vd=vd, label=NN_cluster$label))
  }
  return(scores_whole)
}