source(here::here("methods/outlyingness_scores/Outlyingness_Score.R"))
source(here::here("R/ccds/RK_CCD_New.R"))

# calculate the Outbound Outlyingness Scores of RKCCD
# niter is the number of simulated data set for CSR test if simul is not provided
# d: dimensionality
RKCCD_OOS = function(datax, simul=NULL, d, quant=0.99, tiebreak=getOption("ccd.tiebreak","eq10"), parts=FALSE){
  # The digraph based on RKCCD
  graph = RKCCD_correct_quant(datax,r.seq=10, dom.method="greedy2", 
                                  quan=quant, simul=simul, niter=1000, scores=T)
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

# calculate the Inbound Outlyingness Scores of RKCCD
# niter is the number of simulated data set for CSR test if simul is not provided
# d: dimensionality
RKCCD_IOS = function(datax, simul=NULL, d, quant=0.99, min.cls=0, tiebreak=getOption("ccd.tiebreak","eq10"), parts=FALSE){
  # The digraph based on RKCCD
  RK_cluster = RKCCD_correct_quant(datax,r.seq=10, dom.method="greedy2", 
                              quan=quant, simul=simul, niter=1000, scores=T, min.cls=min.cls)
  # The radius of each point
  R = RK_cluster$R[order(RK_cluster$D)]
  vd = Vic_Den(datax,R,d) # the vicinity density
  # Number of clusters
  n_cls = RK_cluster$si.ind
  member = unique(RK_cluster$label)
  # calculate the scores for each cluster and standardization
  scores = lapply(member,function(x){
    index_cls = which(RK_cluster$label==x)
    cluster = datax[index_cls,]
    R_cls = R[index_cls]
    score_cls = IOS(cluster,R_cls,d)
    score_cls = std_MADN(score_cls)
    return(score_cls)
  })
  scores_whole = rep(0,length(datax[,1]))
  for(i in 1:n_cls){
    scores_whole[which(RK_cluster$label==member[i])] = scores[[i]]
  }
  
  # break ties
  raw_scores = scores_whole
  scores_whole = break_ties(scores_whole, vd, mode = tiebreak)
  if(parts){
    return(list(score=scores_whole, raw=raw_scores, vd=vd, label=RK_cluster$label))
  }
  return(scores_whole)
}