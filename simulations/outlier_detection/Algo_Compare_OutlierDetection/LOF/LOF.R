library(dbscan)
# data: the data set to compute LOF
# L_MinPts and U_MinPts are the upper and lower bound of MinPts
LOF = function(data, L_MinPts=11, U_MinPts=30){
  range_MinPts = c(L_MinPts:U_MinPts)
  
  LOF_mat = sapply(range_MinPts, function(t){
    scores_temp = lof(data, t)
    return(scores_temp)
  })
  
  scores = apply(LOF_mat, 1, function(row){
    return(max(row))
  })
  
  return(scores)
}

# data: the input data set 
# cont: the contamination level
# es_label: the estimated label by an algorithm, 0 represents outlier
# n1 and n0 are the numbers of regular observations and outliers
count_LOF2 = function(n1, n0, es_labels){
  n = n1 + n0
  TPR = sum(es_labels[c((n1+1):n)]==0)/n0
  TNR = sum(es_labels[c(1:n1)]!=0)/n1
  BA = (TNR+TPR)/2
  recall = TPR
  precision = n0*TPR/(n0*TPR+(1-TNR)*(n-n0))
  F2 = 5*precision*recall/(4*precision+recall)
  if(is.na(F2)) F2=0
  return(c(TPR=TPR,TNR=TNR,BA=BA,F2=F2))
}