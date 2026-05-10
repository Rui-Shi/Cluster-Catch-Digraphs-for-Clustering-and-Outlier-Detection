library(isotree)

# Isolation forest for outlier detection
# ntree: the number of trees to construct
# sample_size: the size of each tree
# threshold: threshold for outlier scores
ISO = function(data, ntrees = 1000, sample_size = 256, threshold = 0.55){
  if(sample_size>dim(data)[1]) sample_size=dim(data)[1]
  model = isolation.forest(data, ntrees = ntrees, sample_size = sample_size)
  outlier_scores <- predict(model, data)
  label = ifelse(outlier_scores > threshold, 0, 1)
  return(label)
}

# data: the input data set 
# cont: the contamination level
# es_label: the estimated label by an algorithm, 0 represents outlier
# n1 and n0 are the numbers of regular observations and outliers
count_ISO2 = function(n1, n0, es_labels){
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
