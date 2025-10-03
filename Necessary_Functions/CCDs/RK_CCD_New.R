library(MASS)
library(cluster)
library(igraph)
source("/mmfs1/home/rzs0112/code_working_folder/ccds/ccdfunctions.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/Kest.R")

# ccd clustering that find the dominating set with greedy alg
# translation correction as edge correction
# the boxes are are given with a K function estimation
# a quantile is selected from the simulated envelopes 
# old name = ccd3.clustering_correct_quantile
# datax, the data set
# low.num, lowest number of a box cardinality 
# r.seq, the number of radii used for each window, proximity region
# dom.method, method for finding the dominating set
# quan, quantile to be used for confidence interval
# simul, provided simulations, if null, compute
# method: the way to specify t, non-dynamic or dynamic
rccd.clustering_correct_quantile <- function(datax,low.num,r.seq,method="non-dynamic",
                                             dom.method="greedy2", quan, simul=NULL, niter, scores){
  #low.num=2;r.seq=10;method="non-dynamic";dom.method="greedy2";quan=0.99;niter=100
  # info
  nr <- nrow(datax)
  nc <- ncol(datax)
  ddatax <- as.matrix(dist(datax))
  mac.eps <- .Machine$double.eps
  
  # find all radii, and the dominating set 
  ccd.info <- ccd.Kest.edge.quantile(datax, ddatax, low.num, method, r.seq, quan, simul, niter, scores=scores)
  
  r <- ccd.info$R
  #M <- matrix(as.integer(ddatax < r+mac.eps), length(r))
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
  density[is.infinite(density)]=0
  
  return(list(Int.D=Int.D,Int.R=Int.R,D=D,R=R, catch=catch, density=density))
}


# the validation function for the clusters, does not need the actual cluster labels
# old name = ccd3.clustering.nonvalidddatax <- as.matrix(dist(datax))
# graph is the digraph CCD
# ddatax is the distance matrix
# cls is the actual classes 
rccd.clustering.nonvalid <- function(graph,ddatax){
  
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
# min.cls: the minimum percentage accepted as a cluster
rccd.silhouette <- function(graph, ddatax, cls=NULL, min.cls=0, ind=NULL, lenDlimit=Inf){
  n = nrow(ddatax)
  lenD = length(which(graph$catch>=round(min.cls*n)))
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
      result.temp <- rccd.clustering.nonvalid(dgraph,ddatax)
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
  
  result= as.vector(result)
  return(list(si=maxsi,si.ind=si.ind,label=result))
}

# ccd function that the radii of all points with K function, translation correction as the edge correction, fast version
# incorporates quantiles from the simulated envelopes
# old name = ccd.Kest2.quantile
# low.num is the lowest cardinality of a ball
# r.seq is the number of breaks in the Kest analysis of radii
# quantile to be used as the confidence interval max
# simul is the provided simulation for the problem, if null, compute
# scores: whether to calculate the outlyingness scores
ccd.Kest.edge.quantile <- function(dx, ddx, low.num, method="non-dynamic", r.seq, quan, simul=NULL, niter, scores=F){
  
  n <- nrow(dx)
  d <- ncol(dx)
  R <- rep(0,n)
  #rr <- seq(1/r.seq, 1, 1/r.seq)
  
  # the new simulation function
  if(!is.null(simul)){
    Kest.slopes <- simul
  } else if(method=="non-dynamic"){
    Kest.slopes <- Kest.simpois.edge.quantile(n, d, r.seq, quan, niter)
  } else {
    Kest.slopes <- Kest.simpois.edge.quantile.dynamic(n, d, r.seq, quan, niter)
  }
  rr = Kest.slopes$r
  
  if(!scores){
    for(i in 1:n){
      o.d <- order(ddx[i,]) # the distance vector if i_th object
      for(j in low.num:n){
        r <- ddx[i,o.d[j]]*rr
        sc <- ddx[i,o.d[j]]
        Kest.obs <- Kest.f.edge(ddx[o.d[1:j],o.d[1:j]],r,sc,d)
        
        # check the values, if rejected, set the R[i] as the radius
        lo.hi <- Kest.slopes$quan[[as.character(quan)]][j,]
        flag <- (Kest.obs > lo.hi)
        if(any(flag)){
          if(j==low.num) R[i] <- 0
          else R[i] <- ddx[i,o.d[j-1]]
          break
        } else if(j==n){
          R[i] <- 0
          break
        }
      }
    }
  } else {
    for(i in 1:n){
      o.d <- order(ddx[i,]) # the distance vector if i_th object
      for(j in low.num:n){
        r <- ddx[i,o.d[j]]*rr
        sc <- ddx[i,o.d[j]]
        Kest.obs <- Kest.f.edge(ddx[o.d[1:j],o.d[1:j]],r,sc,d)
        
        # check the values, if rejected, set the R[i] as the radius
        lo.hi <- Kest.slopes$quan[[as.character(quan)]][j,]
        flag <- (Kest.obs > lo.hi)
        if(any(flag)){
          if(j==low.num) R[i] <- 0
          else R[i] <- ddx[i,o.d[j-1]]
          break
        } else if(j==n){
          R[i] <- 0
          break
        }
        if(R[i]==0){R[i]=sort(ddx[i,])[2]} # avoid 0 radius (necessary for outlyingness scores!)
      }
    }
  }
  return(list(R=R,KS=NULL))
}



# Aggregate the functions above into one
RKCCD_correct_quant <- function(datax,low.num=2,r.seq=10,
                                dom.method="greedy2", quan=0.99, simul=NULL, cls=NULL,ind=NULL, lenDlimit=Inf, niter=100, scores=F, min.cls=0){
  
  if(low.num<2) stop("the lowest cardinality of a ball(low.num) can't be less than 2")
  
  if(is.null(lenDlimit)) stop("the maximum index of the dominating point to be checked for silhouette(lenDlimit) can't be NA or NULL")
  
  if(is.na(lenDlimit)) stop("the maximum index of the dominating point to be checked for silhouette(lenDlimit) can't be NA or NULL")
  
  graph <- rccd.clustering_correct_quantile(datax, low.num, r.seq, method="non-dynamic", dom.method, quan, simul, niter,scores=scores)
  ddatax <- as.matrix(dist(datax))
  return(c(graph, rccd.silhouette(graph, ddatax, cls, min.cls=min.cls, ind, lenDlimit)))
}


# ccd clustering for arbitrary shapes of data
# mutual catch graph of RK-CCD
# edge correction included
# the balls are given with a K function estimation
# old name = ccd4.clustering_correct
# datax, the data set
# low.num, lowest number of a box cardinality 
# r.seq, the number of radii used for each window, proximity region
# method, the method for dominating set algorithm
# simul, the simulation results if provided
rccd.clustering.mutual.connected_correct <- function(datax,low.num,r.seq, method="non-dynamic",
                                                     dom.method=NULL, quan, simul=NULL, niter,scores=F){
  # info
  nr <- nrow(datax)
  nc <- ncol(datax)
  ddatax <- as.matrix(dist(datax))
  mac.eps <- .Machine$double.eps
  
  # find all radii, and the dominating set 
  ccd.info <- ccd.Kest.edge.quantile(datax, ddatax, low.num, method, r.seq, quan, simul, niter, scores=scores)
  
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
  
  # find the edges of the mutual catch graph of RK-CCDs
  M <- M[D,D]
  M.dom <- diag(T,nrow(M))
  for(i in 1:(nrow(M.dom)-1)){
    for(j in (i+1):nrow(M.dom)){
      temp <- (M[i,j] & M[j,i])
      if(any(temp)) M.dom[i,j] <- M.dom[j,i] <- T
    }
  }
  
  # find the disconnected components that are the clusters
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
rccd.clustering.nonvalid.knn <- function(Mgraph,datax,k = NULL,cls=NULL){
  # ddatax = Mgraph$dist
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


# for UN-MCCD
# the validation function for the clusters, does not need the actual cluster labels
# k is the number of nearest neighbors to consider, if not specified, k = min(dimension, data size)
# graph is the mutual catch graph of RK-CCD
# datax is the data set
# ddatax is the distance matrix
# cls: number of connected cluster to consider
rccd.clustering.nonvalid.knn1 <- function(graph,datax,k = NULL,cls=NULL){
  ddatax = graph$dist.M
  # datax = datax[graph$D,]
  # ddatax = graph$dist[graph$D,graph$D]
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
rccd.silhouette_mutual <- function(Mgraph,datax,ind=NULL, lenClimit=Inf, k=NULL, min.cls = 0){
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
      label = rccd.clustering.nonvalid.knn(Mgraph.temp,datax,k,cls=NULL)
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


# For UN-MCCD
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
rccd.silhouette_mutual1 <- function(graph, datax,ind=NULL, lenClimit=Inf, k=NULL, min.cls = 0){
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
      label = rccd.clustering.nonvalid.knn1(graph.temp,datax,k,cls=NULL)
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