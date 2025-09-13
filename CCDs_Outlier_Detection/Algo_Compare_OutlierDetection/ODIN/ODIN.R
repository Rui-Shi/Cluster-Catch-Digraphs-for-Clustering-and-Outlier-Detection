library(FNN)

# Outlier Detection Using Indegree Number (ODIN)

ODIN <- function(data, k=NULL, indegree_threshold=NULL, quant=NULL) {
  # Calculate kNN graph
  if(is.null(k)){k=round(sqrt(dim(data)[1]))}
  if(is.null(indegree_threshold)){indegree_threshold=round(dim(data)[1]^(1/3))}
  knn_graph <- get.knn(data, k = k)
  
  # Calculate indegree (Corrected Calculation)
  indegree <- sapply(1:nrow(data), function(i) sum(i == knn_graph$nn.index))
  
  if(!is.null(quant)){indegree_threshold = quantile(indegree, quant)}
  # Identify outliers based on threshold
  outliers <- indegree <= indegree_threshold
  outliers <- ifelse(outliers, 0, 1)
  return(outliers)
}

# data: the input data set 
# n: the size of entire data set
# cont: the contamination level
# es_label: the estimated label by an algorithm, 0 represents outlier
count_ODIN = function(n, cont, es_labels){
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
count_ODIN2 = function(n1, n0, es_labels){
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