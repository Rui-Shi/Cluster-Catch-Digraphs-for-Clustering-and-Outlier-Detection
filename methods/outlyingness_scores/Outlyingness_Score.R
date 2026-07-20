# Calculate the vicinity density of each points
# data: the entire data set or cluster
# the radius of each point
# d: dimensionality
Vic_Den = function(data, R, d){
  n = dim(data)[1]
  ddatax = as.matrix(dist(data))
  if(is.null(n)) {n=1; ddatax=matrix(0,nrow=1,ncol=1)} # the case when there is only one observation
  Den = sapply(1:n, function(x){
    size = length(which(ddatax[x,]<=R[x]))
    # return(size/(R[x]^d))
    # return((size/R[x]^d)^(1/d))  # overflows: R^d = Inf for R > exp(709/d)
    #                              # (~13.3 at d=274), giving Den = 0 -> OOS NaN,
    #                              # IOS Inf. See revision_experiments/
    #                              # test_outlyingness_density.R (2026-07-16).
    return(size^(1/d)/R[x])       # algebraically identical, numerically stable
  })
  return(Den)
}

# Calculate the Outbound Outlying-ness Score (OOS) of each points
# data: the entire data set or cluster
# the radius of each point
# d: dimensionality
OOS = function(data, R, d){
  n = dim(data)[1]
  ddatax = as.matrix(dist(data))
  if(is.null(n)) {n=1; ddatax=matrix(0,nrow=1,ncol=1)} # the case when there is only one observation
  diag(ddatax) = Inf
  Den = Vic_Den(data, R, d)
  scores = sapply(1:n,function(x){
    out_nei.index = which(ddatax[x,]<=R[x])
    score = mean(Den[out_nei.index])/Den[x]
  })
  return(scores)
}

# Calculate the Inbound Outlying-ness Score (OOS) of each points
# data: the entire data set or cluster
# the radius of each point
# d: dimensionality
IOS = function(data, R, d){
  n = dim(data)[1]
  ddatax = as.matrix(dist(data))
  if(is.null(n)) {n=1; ddatax=matrix(0,nrow=1,ncol=1)} # the case when there is only one observation
  Den = Vic_Den(data, R, d)
  scores = sapply(1:n,function(x){
    in_nei.index = which(ddatax[x,]<=R)
    score = 1/sum(Den[in_nei.index])
  })
  return(scores)
}

# standardization using MADN, with the manuscript's degenerate-case fallback
# chain MADN -> SD -> 0 (CCDwScores.tex line 625, item B3): when MADN = 0 the
# denominator is replaced by the sample SD (median centering unchanged); if SD
# also vanishes (or is undefined, i.e. a singleton cluster), all scores are 0.
# See revision_experiments/test_std_madn.R.
# x: a set of sample point
std_MADN = function(x){
  if(mad(x)!=0){
    # s = abs((x-median(x))/mad(x))
    s = (x-median(x))/mad(x) # allow negative scores
  } else {
    sdx = sd(x)
    if(!is.na(sdx) && sdx!=0){
      s = (x-median(x))/sdx  # MADN=0 fallback: replace denominator with SD
    } else {
      s = rep(0,length(x))   # MADN=SD=0 (or singleton): degenerate cluster
    }
  }
  return(s)
}