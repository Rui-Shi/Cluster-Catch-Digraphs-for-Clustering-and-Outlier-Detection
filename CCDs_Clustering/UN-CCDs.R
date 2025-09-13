source("/media/rui/exNVME/code_working_folder/ccds/UN_CCD.R")

# main function for outlier function
# datax: the data set
# simul: the simulated windows or envelope 
# cls: the actual classes 
# ind: the set of indices of clusterings to be checked
# min.cls: the minimum percentage accepted as a cluster

UNCCD_clustering = function(datax, simul, method="ascend", niter=1000, cls=NULL, ind=NULL, min.cls=0){
  UN_graph <- nnccd_clustering_quantile(datax, low.num=3, quantile="lower", 
                                     method=method, dom.method="greedy2", 
                                     simul=simul, niter=niter)
  UN_result = nnccd.silhouette(UN_graph,datax,cls=NULL,ind=NULL, lenDlimit=Inf, min.cls=min.cls)
  UN_result = c(UN_graph,UN_result)
  
  # a list of clusters of the dataset
  UN_result$cluster = lapply(1:UN_result$si.ind, function(x){return(datax[UN_result$label==x,])})
  
  return(UN_result)
}