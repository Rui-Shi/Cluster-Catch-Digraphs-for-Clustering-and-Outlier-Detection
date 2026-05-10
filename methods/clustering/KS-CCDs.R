source("/media/rui/exNVME/code_working_folder/ccds/ccd_ks_NEW.R")
source("/media/rui/exNVME/code_working_folder/ccds/ccdfunctions.R")

# datax: the input data set
# m: K-S statistics parameter
# cls: the actual classes 
# ind: the set of indices of clusterings to be checked
# lenDlimit: the maximum index of the dominating point to be checked for silhouette
# min.cls: the minimum percentage accepted as a cluster
# low.num: the min number of each covering ball
KSCCD_clustering = function(datax, m=10, cls=NULL, ind=NULL){
  KS_result = ksccd.clustering(datax=datax, m=m, sequential=FALSE,
                               dom.method="greedy2", alpha=0.05)
  
  ddatax <- as.matrix(dist(datax))
  clustering_result = ksccd.silhouette(KS_result, ddatax, cls=NULL, min.cls=0, ind=NULL, lenDlimit=Inf)
  
  KS_result = c(KS_result, clustering_result)
  KS_result$clusters = lapply(1:KS_result$si.ind, function(x){return(datax[KS_result$label==x,])})
  
  return(c(KS_result, clustering_result))
}