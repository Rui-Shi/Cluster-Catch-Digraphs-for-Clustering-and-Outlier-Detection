# source("/mmfs1/home/rzs0112/code_working_folder/ccds/RK_CCD_New.R")
source("/media/rui/exNVME/code_working_folder/ccds/RK_CCD_New.R")

# datax: the input data set
# Simul: the upper and lower envelope of the K-function
# niter: the number of simulations to get Simul
# cls: the actual classes 
# ind: the set of indices of clusterings to be checked
# lenDlimit: the maximum index of the dominating point to be checked for silhouette
# min.cls: the minimum percentage accepted as a cluster
RKCCD_clustering = function(datax, simul=NULL, quant=0.99, cls=NULL,ind=NULL, niter=1000){
  RK_result = RKCCD_correct_quant(datax,r.seq=10, dom.method="greedy2", 
                                  quan=quant, simul=simul, cls=NULL,ind=NULL, niter=niter)
  # a list of clusters of the dataset
  RK_result$clusters = lapply(1:RK_result$si.ind, function(x){return(datax[RK_result$label==x,])})
  
  return(RK_result)
}