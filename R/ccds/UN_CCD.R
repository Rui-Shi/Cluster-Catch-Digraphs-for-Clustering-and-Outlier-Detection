library(MASS)
library(cluster)
library(igraph)
source(here::here("R/ccds/ccdfunctions.R"))
source(here::here("R/ccds/NN_Dist_Est.R"))

# ccd clustering that find the dominating set with greedy alg
# old name = ccd3.clustering_correct_quantile
# datax, the data set
# low.num, lowest number of a box cardinality 
# r.seq, the number of radii used for each window, proximity region
# dom.method, method for finding the dominating set
# quan, quantile to be used for confidence interval
# simul, provided simulations
nnccd_clustering_quantile <- function(datax,low.num=3, quantile="lower", 
                                      method="ascend", dom.method="greedy2", quant=0.90, simul=NULL, niter=100, scores=F){
  # info
  nr <- nrow(datax)
  nc <- ncol(datax)
  ddatax <- as.matrix(dist(datax))
  mac.eps <- .Machine$double.eps
  
  # find all radii, and the dominating set 
  ccd.info <- nnccd.radi(dx=datax, quantile, method, low.num, quant, simul=simul, niter, scores=scores)
  
  r <- ccd.info$R
  # M <- matrix(as.integer(ddatax < r+mac.eps), length(r))
  M <- matrix(as.integer(ddatax <= r), length(r))
  score <- rowSums(M)
  
  # domination method is provided, find the dominating set
  # if not, take all radii
  if(!is.null(dom.method)){
    if(dom.method=="greedy") D <- dominate.mat.greedy(M)
    if(dom.method=="greedy2") D <- dominate.mat.greedy2(M)
    if(dom.method=="ks") D <- dominate.mat.ks(M,ks)
  } else {
    D <- 1:nr
  }
  R <- r[D]
  
  # find the intersection graph and its dominating set
  # the intersection graph is based on whether two radii catched the same points
  MD <- M[D,]
  M.dom <- diag(T,nrow(MD))
  for(i in 1:(nrow(M.dom)-1)){
    for(j in (i+1):nrow(M.dom)){
      temp <- (MD[i,] & MD[j,])
      if(any(temp)) M.dom[i,j] <- M.dom[j,i] <- T
    }
  }
  
  # dominating set of the intersection graph
  D.ind <- dominate.mat.ks(M.dom,score[D])
  
  # keep the 2nd layer of dominating sets 
  Int.D=D[D.ind]
  Int.R=R[D.ind]
  
  # get the catching number and density of all Dominated graphs 
  MDInt <- M[Int.D,,drop=FALSE] #with drop=FALSE, it output a matrix or dataframe rather than a vector
  # catch <- apply(MDInt,1,function(x){
  #               temp <- ddatax[as.logical(x),as.logical(x),drop = FALSE]
  #               diag(temp) <- 0
  #               return(sum(temp)/(sum(x)^2-sum(x)))
  # })
  catch <- rowSums(MDInt)
  density <- catch/(Int.R^nc)
  
  return(list(Int.D=Int.D,Int.R=Int.R,D=D,R=R, catch=catch, density=density))
}


nnccd.clustering.nonvalid <- function(graph,ddatax){
  
  D <- graph$Int.D
  R <- graph$Int.R
  cl <- 1:length(D)
  #print(cl)
  ddx <- matrix(ddatax[,D],ncol=length(D))
  
  # which point is in which ball of dominators
  ddx <- ddx + .Machine$double.eps
  result <- t(apply(ddx,1,function(x){
    ind <- which.min(x/R)
    return(cl[ind])
  }))
  return(result)
}


# incrementally add cluster points to find the clustering with maximum silhouette
# old name ccd3.silhouette
# graph is the digraph CCD
# ddatax is the distance matrix
# cls is the actual classes 
# ind is the set of indices of clusterings to be checked
# lenDlimit is the maximum index of the dominating point to be checked for silhouette
nnccd.silhouette <- function(graph, datax, cls=NULL, min.cls=0, ind=NULL, lenDlimit=Inf){
  ddatax = as.matrix(dist(datax))
  n = nrow(ddatax)
  lenD = length(which(graph$catch>=round(min.cls*n)))
  #lenD <- length(graph$Int.D)
  if(lenDlimit < lenD) lenD <- lenDlimit
  dgraph <- graph
  
  if(is.null(ind)) ind <- 1:nrow(ddatax)
  
  maxsi <- 0
  si.ind <- 1
  result <- rep(1,nrow(ddatax))
  if(lenD<2){
    si.ind=1;result = rep(1,n);maxsil=NULL
  } else {
    for(i in 2:lenD){
      dgraph$Int.D <- graph$Int.D[1:i]
      dgraph$Int.R <- graph$Int.R[1:i]
      #result <- ccd3.clustering.nonvalid(dgraph,ddatax)
      result.temp <- nnccd.clustering.nonvalid(dgraph,ddatax)
      if(length(unique(result.temp[ind])) < 2) {
        datasi <- 0
      } else {
        datasi <- silhouette(result.temp[ind],ddatax[ind,ind])
        datasi <- mean(datasi[,3])  
        #print(datasi)
      }
      if(datasi > maxsi & sum(result.temp==i)>1){
        result <- result.temp
        maxsi <- datasi
        si.ind <- i
      } 
      # print(paste("sil", i," = ", datasi,sep=""))
    }
  }
  result = as.vector(result)
  return(list(si=maxsi,si.ind=si.ind,label=result))
}


# incrementally add cluster points to find the clustering with maximum silhouette
# graph is the mutual catch graph of NN-CCD
# ddatax is the distance matrix
# cls is the actual classes 
# ind is the set of indices to be checked for clusters
# lenClimit is the maximum index of cluster to be checked for silhouette
# ind.C: optimal number of clusters
# label: clustering result
# max.sil: maximum silhouette index obtained
# min.cls: the minimum percentage accepted as a cluster
nnccd.silhouette_mutual <- function(Mgraph,datax,ind=NULL, lenClimit=Inf, k=NULL, min.cls = 0){
  n = nrow(datax)
  cls.size = sort(table(Mgraph$D.member),decreasing = T)
  
  # ddatax = Mgraph$dist
  ddatax = Mgraph$dist[Mgraph$D,Mgraph$D]
  lenC = min(length(Mgraph$member),lenClimit)
  if(is.null(ind)){ind = n} else {ind = min(ind, n)}
  
  lenC = length(which(cls.size>=round(min.cls*n)))
  
  if(lenC<2){
    sil.ind=1;label = rep(Mgraph$member[1],n);maxsil=NULL
  } else {
    result = lapply(2:lenC,function(t){
      Mgraph.temp = Mgraph
      Mgraph.temp$member = Mgraph$member[1:t]
      label = nnccd.clustering.nonvalid.knn(Mgraph.temp,datax,k,cls=NULL)
      sil = mean(silhouette(label,ddatax)[1:ind,3])
      return(list(label=label,sil=sil))
    })
    
    sil.ind = NULL
    maxsil = -1e20
    label = NULL
    for(i in 1:length(result)){
      if(result[[i]]$sil>maxsil){sil.ind = i+1;label = result[[i]]$label;maxsil = result[[i]]$sil}
    }
  }
  return(list(sil=maxsil,sil.ind=sil.ind,label=label))
}

# For SUN-MCCD
# incrementally add cluster points to find the clustering with maximum silhouette
# graph is the mutual catch graph of RK-CCD
# ddatax is the distance matrix
# cls is the actual classes 
# ind is the set of indices of clusterings to be checked
# lenClimit is the maximum index of cluster to be checked for silhouette
# ind.C: optimal number of clusters
# label: clustering result
# max.sil: maximum silhouette index obtained
# min.cls: the minimum percentage accepted as a cluster
nnccd.silhouette_mutual1 <- function(graph, datax,ind=NULL, lenClimit=Inf, k=NULL, min.cls = 0){
  n = nrow(datax)
  cls.size = sort(table(graph$member),decreasing = T)
  
  ddatax = graph$dist.M
  lenC = min(length(graph$D.member),lenClimit)
  if(is.null(ind)){ind = n} else {ind = min(ind, n)}
  
  lenC = length(which(cls.size>=round(min.cls*n)))
  
  if(lenC<2){
    sil.ind=1;label = rep(graph$D.member[1],n);maxsil=NULL
  } else {
    result = lapply(2:lenC,function(t){
      graph.temp = graph
      graph.temp$D.member = graph$D.member[1:t]
      label = nnccd.clustering.nonvalid.knn1(graph.temp,datax,k,cls=NULL)
      sil = mean(silhouette(label,ddatax)[1:ind,3])
      return(list(label=label,sil=sil))
    })
    
    sil.ind = NULL
    maxsil = -1e20
    label = NULL
    for(i in 1:length(result)){
      if(result[[i]]$sil>maxsil){sil.ind = i+1;label = result[[i]]$label;maxsil = result[[i]]$sil}
    }
  }
  return(list(sil=maxsil,sil.ind=sil.ind,label=label))
}


# dx: data set
# method: "upper" : upper tail, "lower": lower tail, "two": two sided tails
# ccd function that compute the radii of all points with K function
# incorporates quantiles from the simulated envelopes
# low.num is the lowest cardinality of a ball
# quantile to be used as the confidence interval max
# simul is the provided simulation for the problem, if null, compute
# scores: whether to calculate the outlyingness scores
# Incremental nearest-neighbour update (2026-09-23).
# The earlier version rebuilt every nearest-neighbour distance of the candidate
# ball from a fresh submatrix at each candidate radius (NNDest.dist.f), which
# costs O(j^2) at a ball of j points and O(n^4) over the whole search in the
# worst case. Here the within-ball NN distances are kept in a vector `nn` in the
# same member order as before and updated when one point enters (ascending) or
# leaves (descending) the ball, so each step costs O(j). mean() and median() are
# still taken over the same values in the same order, so the statistics, and
# therefore the radii, are identical to the earlier version
# (checked by revision_experiments/tr1/92b_validate_incremental_radi.R against
# R/ccds/UN_CCD_radi_recompute_reference.R).

# Ascending search for centre i: ball members o.d[2:j]; returns the radius.
nnccd.radi.ascend.one <- function(i, ddx, o.d, n, low.num, env.ave, env.med) {
  m0 <- o.d[2:low.num]                       # members at the first step
  if (length(m0) < 2) {                      # degenerate ball: defer to the reference statistic
    for (j in low.num:n) {
      r <- ddx[i, o.d[j]]
      s <- NNDest.dist.f(ddx[o.d[2:j], o.d[2:j]], r)
      if (s$averge < env.ave[j - 1] | s$median < env.med[j - 1]) {
        return(if (j == low.num) 0 else ddx[i, o.d[j - 1]])
      }
    }
    return(0)
  }
  sub <- ddx[m0, m0, drop = FALSE]; diag(sub) <- Inf
  nn <- apply(sub, 1, min)
  for (j in low.num:n) {
    if (j > low.num) {                        # point o.d[j] enters the ball
      p <- o.d[j]
      dp <- ddx[p, o.d[2:(j - 1)]]
      nn <- c(pmin(nn, dp), min(dp))
    }
    r <- ddx[i, o.d[j]]
    if (mean(nn) / r < env.ave[j - 1] | median(nn) / r < env.med[j - 1]) {
      return(if (j == low.num) 0 else ddx[i, o.d[j - 1]])
    }
  }
  0
}

# Descending search for centre i: ball members o.d[j:(n-1)]; returns the radius.
# nn1.val/nn1.idx/nn2.val: each point's nearest and second-nearest distance to
# any other point (computed once per data set), used to start the full ball in O(n).
nnccd.radi.descend.one <- function(i, ddx, o.d, n, low.num, env.ave.rev, env.med.rev,
                                   nn1.val, nn1.idx, nn2.val) {
  mem <- o.d[1:(n - 1)]                       # members at j = 1
  e <- o.d[n]                                 # the one point outside the ball
  nn <- ifelse(nn1.idx[mem] == e, nn2.val[mem], nn1.val[mem])   # position k <-> point mem[k]
  for (j in 1:(n - low.num)) {
    if (j > 1) {                              # point o.d[j-1] leaves; members are mem[j:(n-1)]
      x <- o.d[j - 1]
      pos <- j:(n - 1)
      hit <- pos[nn[pos] == ddx[mem[pos], x]]
      for (k in hit) {
        others <- mem[pos[pos != k]]
        nn[k] <- if (length(others)) min(ddx[mem[k], others]) else Inf
      }
    }
    r <- ddx[i, o.d[j]]
    cur <- nn[j:(n - 1)]
    if (mean(cur) / r > env.ave.rev[j + 2] & median(cur) / r > env.med.rev[j + 2]) {
      return(r)
    }
  }
  0
}

nnccd.radi <- function(dx, quantile="lower", method="ascend", low.num, quant, simul=NULL, niter, scores=F){

  ddx <- as.matrix(dist(dx)) # the distance matrix
  n <- nrow(dx)
  d <- ncol(dx)
  R <- rep(0,n)

  if(quantile=="lower"){
    if(!is.null(simul)) {NN.envelop <- list(average=simul$average[1:n],median=simul$median[1:n])}
    else {NN.envelop <- NNDest.simpois.lower.quant(n, d, quant, niter)}
    if(method=="descend"){
      env.ave.rev <- rev(NN.envelop$average); env.med.rev <- rev(NN.envelop$median)
      dd0 <- ddx; diag(dd0) <- Inf
      nn1.idx <- apply(dd0, 1, which.min)
      nn1.val <- dd0[cbind(1:n, nn1.idx)]
      dd0[cbind(1:n, nn1.idx)] <- Inf
      nn2.val <- apply(dd0, 1, min)
      rm(dd0)
    }
    for(i in 1:n){
      if(method == "ascend"){
        o.d <- order(ddx[i,]) # the ascending distance order for i_th object
        R[i] <- nnccd.radi.ascend.one(i, ddx, o.d, n, low.num, NN.envelop$average, NN.envelop$median)
      }
      if(method=="descend"){
        o.d <- order(ddx[i,], decreasing=T) # the descending distance order for i_th object
        R[i] <- nnccd.radi.descend.one(i, ddx, o.d, n, low.num, env.ave.rev, env.med.rev,
                                       nn1.val, nn1.idx, nn2.val)
      }
      if(scores && R[i]==0){R[i]=sort(ddx[i,])[2]} # avoid 0 radius (necessary for outlyingness scores!)
    }
  }
  return(list(R=R,KS=NULL))
}


# nnccd clustering for arbitrary shapes of data
# mutual catch graph of nn-CCD
# datax, the data set
# low.num, lowest number of a box cardinality 
# quant: the critial quantile for the CSR test based on NN-CCD
# method, the method for dominating set algorithm
# simul, the simulation results if provided
nnccd.clustering.mutual.connected <- function(datax,low.num=3, quantile="lower", 
                                              method="ascend", dom.method="greedy2", quant=0.90, simul=NULL, niter=100, scores=F){
  # info
  nr <- nrow(datax)
  nc <- ncol(datax)
  ddatax <- as.matrix(dist(datax))
  mac.eps <- .Machine$double.eps
  
  # find all radii, and the dominating set 
  ccd.info <- nnccd.radi(dx=datax, quantile, method, low.num, quant, simul=simul, niter, scores=scores)
  
  r <- ccd.info$R
  # M <- matrix(as.integer(ddatax < r + mac.eps), length(r))
  M <- matrix(as.integer(ddatax <= r), length(r))
  diag(M) <- 1
  score <- rowSums(M)
  
  # domination method is provided, find the dominating set
  # if not, take all radii
  if(!is.null(dom.method)){
    if(dom.method=="greedy") D <- dominate.mat.greedy(M)
    if(dom.method=="greedy2") D <- dominate.mat.greedy2(M)
    if(dom.method=="ks") D <- dominate.mat.ks(M,ks)
  } else {
    D <- 1:nr
  }
  R <- r[D]
  
  # find the edges of the mutual catch graph of NN-CCDs
  M <- M[D,D]
  M.dom <- diag(T,nrow(M))
  for(i in 1:(nrow(M.dom)-1)){
    for(j in (i+1):nrow(M.dom)){
      temp <- (M[i,j] & M[j,i])
      if(any(temp)) M.dom[i,j] <- M.dom[j,i] <- T
    }
  }
  
  # find the connected components that are the clusters
  D.member <- components.mat(M.dom)
  member <- order(table(D.member),decreasing = T)
  return(list(member = member, D.member=D.member, D=D, R=R, dist = ddatax))
}


# calculate the local distance-based outlier factor(LDOF) of a point x to a cluster
# x: a given point, cluster: a given cluster of points
# k is the number of nearest neighbors to consider, if not specified, k = min(dimension, data size)
# dist[[1]] is the dist vector from the point to the cluster, dist[[2]] is the dist matrix of the cluster
LDOF <- function(x, cluster, k = NULL, dist){
  
  if(is.null(k)) {
    k = min(ncol(cluster), nrow(cluster))
  } else {
    k = min(nrow(cluster), k)
  }
  
  ind = order(dist[[1]])[1:k]
  ave1 =  mean(dist[[1]][ind])
  
  if(k==1){ave2 = 1e-20} else {
    inner.dist = dist[[2]][ind,ind]
    ave2 = mean(inner.dist[which(upper.tri(inner.dist))])
  }
  
  result = ave1/ave2
  return(result)
}


# the validation function for the clusters, does not need the actual cluster labels
# k is the number of nearest neighbors to consider, if not specified, k = min(dimension, data size)
# graph is the mutual catch graph of RK-CCD
# datax is the data set
# ddatax is the distance matrix
# cls: number of connected cluster to consider
nnccd.clustering.nonvalid.knn <- function(Mgraph,datax,k = NULL,cls=NULL){
  #ddatax = Mgraph$dist
  datax = datax[Mgraph$D,]
  ddatax = Mgraph$dist[Mgraph$D,Mgraph$D]
  if(is.null(cls)){cls = length(Mgraph$member)}
  n = nrow(datax)
  d = ncol(datax)
  cluster.ind = lapply(1:cls,function(x){
    ind = which(Mgraph$D.member==Mgraph$member[x])
    return(ind)
  })
  
  # which point is in which cluster
  result = sapply(1:n, function(a){
    if(!any(Mgraph$D.member[a]==Mgraph$member)){
      dist.list = sapply(1:cls, function(b){
        dist = list(ddatax[a,cluster.ind[[b]]],ddatax[cluster.ind[[b]],cluster.ind[[b]]])
        x = datax[a,]
        cluster = matrix(datax[cluster.ind[[b]],],ncol = d)
        dist.cls = LDOF(x, cluster, k, dist)
        return(dist.cls)
      })
      ind.cls = which(dist.list==min(dist.list))
      return(Mgraph$member[ind.cls])
      
    } else {return(Mgraph$D.member[a])}
  })
  return(result)
}


# For SUN-MCCD
# the validation function for the clusters, does not need the actual cluster labels
# k is the number of nearest neighbors to consider, if not specified, k = min(dimension, data size)
# graph is the mutual catch graph of RK-CCD
# datax is the data set
# ddatax is the distance matrix
# cls: number of connected cluster to consider
nnccd.clustering.nonvalid.knn1 <- function(graph,datax,k = NULL,cls=NULL){
  ddatax = graph$dist.M
  if(is.null(cls)){cls = length(graph$D.member)}
  n = nrow(datax)
  d = ncol(datax)
  cluster.ind = lapply(1:cls,function(x){
    ind = which(graph$member==graph$D.member[x])
    return(ind)
  })
  
  # which point is in which cluster
  result = sapply(1:n, function(a){
    if(!any(graph$member[a]==graph$D.member)){
      dist.list = sapply(1:cls, function(b){
        dist = list(ddatax[a,cluster.ind[[b]]],ddatax[cluster.ind[[b]],cluster.ind[[b]]])
        x = datax[a,]
        cluster = matrix(datax[cluster.ind[[b]],],ncol = d)
        dist.cls = LDOF(x, cluster, k, dist)
        return(dist.cls)
      })
      ind.cls = which(dist.list==min(dist.list))
      return(graph$D.member[ind.cls])
      
    } else {return(graph$member[a])}
  })
  return(result)
}
