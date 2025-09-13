# count the success rate and false positive rate
count = function(x){
  clsNum = length(outlier.result[[x]]$label)
  outlier.true = matrix(data.list[[x]][-c(1:(n-n0)),], ncol=length(data.list[[x]][1,])) # true outliers
  outlier_detect = NULL # outliers captured
  for(i in 1:clsNum){
    mode = which.max(table(outlier.result[[x]]$label[[i]]))
    outlier_temp = matrix(outlier.result[[x]]$clusters[[i]],ncol=d)[outlier.result[[x]]$label[[i]]!=mode,]
    outlier_detect = rbind(outlier_detect,outlier_temp)
  }
  success_rate = length(intersect(outlier.true[,1], outlier_detect[,1]))/n0
  false_positive = max((nrow(outlier_detect)-success_rate*n0)/(n-n0),0)
  recall = success_rate
  precision = n0*success_rate/(n0*success_rate+false_positive*(n-n0))
  F2 = 5*precision*recall/(4*precision+recall)
  BA=(success_rate+1-false_positive)/2
  if(is.na(F2)) F2=0
  return(c(success_rate=success_rate, true_positive=1-false_positive, BA=BA, F2_score = F2))
}

# count the success rate and false positive rate for clx_cls.R
count1 = function(x){
  n = sum(data.num[[x]])
  n0 = data.num[[x]][2]
  clsNum = length(outlier.result[[x]]$label)
  outlier.true = matrix(data.list[[x]][-c(1:(n-n0)),], ncol=length(data.list[[x]][1,])) # true outliers
  outlier_detect = NULL # outliers captured
  for(i in 1:clsNum){
    mode = which.max(table(outlier.result[[x]]$label[[i]]))
    outlier_temp = matrix(outlier.result[[x]]$clusters[[i]],ncol=d)[outlier.result[[x]]$label[[i]]!=mode,]
    outlier_detect = rbind(outlier_detect,outlier_temp)
  }
  success_rate = length(intersect(outlier.true[,1], outlier_detect[,1]))/n0
  false_positive = max((nrow(outlier_detect)-success_rate*n0)/(n-n0),0)
  recall = success_rate
  precision = n0*success_rate/(n0*success_rate+false_positive*(n-n0))
  F2 = 5*precision*recall/(4*precision+recall)
  if(is.na(F2)) F2=0
  return(c(success_rate=success_rate,false_positive=false_positive,F2_score = F2))
}


# count FP and SR of Heuristic 1
count2 = function(x){
  outlier.index = c(1:n)[-c(1:(n-n0))] # the index of outliers
  usual.index = c(1:n)[c(1:(n-n0))] # the index of usual data points\
  members = components_result[[x]]$member
  
  mode = which.max(table(members))# the majority labels
  FR = length(which(members[usual.index]!=mode))/n1
  SR = length(which(members[outlier.index]!=mode))/n0
  return(c(success_rate=SR,false_positive=FR))
}


# count TNR and TPR of outlyingness score
# threshold: the threshold for outlyingness scores.
# n: number of observations
# n0: number of outliers
# scores: outlyingness scores
# x: index
count_scores = function(x, scores, threshold, n, n0){
  label = rep(0,n)
  score = scores[[x]]
  label[which(score>=threshold)]=2 # outliers are labeled as 2
  label[which(score<threshold)]=1 # regular observations are labeled as 1
  TNR = length(which(label[1:(n-n0)]==1))/(n-n0)
  TPR = length(which(label[(n-n0+1):n]==2))/n0
  BA = (TNR+TPR)/2
  recall = TPR
  precision = n0*TPR/(n0*TPR+(1-TNR)*(n-n0))
  F2 = 5*precision*recall/(4*precision+recall)
  if(is.na(F2)) F2=0
  return(c(TPR=TPR,TNR=TNR,BA=BA,F2=F2))
}

count_scores1 = function(x, scores, threshold){
  n = sum(data.num[[x]])
  n0 = data.num[[x]][2]
  label = rep(0,n)
  score = scores[[x]]
  label[which(score>=threshold)]=0 # outliers are labeled as 0
  label[which(score<threshold)]=1 # regular observations are labeled as 1
  TNR = length(which(label[1:(n-n0)]==1))/(n-n0)
  TPR = length(which(label[(n-n0+1):n]==0))/n0
  BA = (TNR+TPR)/2
  recall = TPR
  precision = n0*TPR/(n0*TPR+(1-TNR)*(n-n0))
  F2 = 5*precision*recall/(4*precision+recall)
  if(is.na(F2)) F2=0
  return(c(TPR,TNR,BA,F2))
}


# Y is the true label (1 is the regular point, 0 is outlier)
# score are the outlier scores
# threshold is the threshold for outlier
count_scores2 = function(Y, score, threshold){
  n = length(Y)
  n0 = sum(Y==0)
  label_pred = rep(0,n)
  label_pred[which(score>=threshold)]=0 # outliers are labeled as 0
  label_pred[which(score<threshold)]=1 # regular observations are labeled as 1
  TNR = length(which(label_pred[1:(n-n0)]==1))/(n-n0)
  TPR = length(which(label_pred[(n-n0+1):n]==0))/n0
  BA = (TNR+TPR)/2
  recall = TPR
  precision = n0*TPR/(n0*TPR+(1-TNR)*(n-n0))
  F2 = 5*precision*recall/(4*precision+recall)
  if(is.na(F2)) F2=0
  return(c(TPR,TNR,BA,F2))
}