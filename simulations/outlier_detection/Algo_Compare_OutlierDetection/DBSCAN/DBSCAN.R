# data: the input data set
# k:  The minimum number of points required within the eps radius to form a dense region (MinPts)
# quant: the (expected) percentage of outliers
library(dbscan)
DBSCAN = function(data, k, quant){
  dist_M = as.matrix(dist(data))
  
  k_dist = apply(dist_M, 1, function(row){  # compute the k-distance for each point
    return(sort(row)[k+1])
  })
  
  eps = quantile(k_dist, 1-quant) # find the eps of DBSCAN based on the outlier level
  
  label = dbscan(data, eps, k)$cluster # conduct DBSCAN
  return(labels=label)
}

# data: the input data set 
# n: the size of entire data set
# cont: the contamination level
# es_label: the estimated label by an algorithm, 0 represents outlier
count_DBSCAN = function(n, cont, es_labels){
  n0 = round(n*cont)
  n1 = n-n0
  TPR = sum(es_labels[c((n1+1):n)]==0)/n0
  TNR = sum(es_labels[c(1:n1)]!=0)/n1
  BA = (TNR+TPR)/2
  recall = TPR
  precision = n0*TPR/(n0*TPR+(1-TNR)*(n-n0))
  F2 = 5*precision*recall/(4*precision+recall)
  if(is.na(F2)) F2=0
  return(c(TPR=TPR,TNR=TNR,BA=BA,F2=F2))
}


# data: the input data set 
# cont: the contamination level
# es_label: the estimated label by an algorithm, 0 represents outlier
# n1 and n0 are the numbers of regular observations and outliers
count_DBSCAN2 = function(n1, n0, es_labels){
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