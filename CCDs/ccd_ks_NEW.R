# Codes of KS-CCD clustering
# ccd clustering that find the dominating set with greedy algorithm
# old name = ccd1.clustering
# datax, the data set
# m, K-S statistics parameter
# sequential, logical whether to use sequential K-S or not
# dom.method: 
#     greedy, dominating sets greedy on outdegree
#     ks    , ominating sets greedy on ks statistics
ksccd.clustering <- function(datax,m,sequential=FALSE,
                             dom.method="greedy2",alpha=0.05){
  
  # info
  nr <- nrow(datax)
  nc <- ncol(datax)
  ddatax <- as.matrix(dist(datax))
  # ddatax <- ddatax/max(ddatax)
  # mac.eps <- .Machine$double.eps ### A number that is extremely small ### 
  
  # find all radii, and the dominating set 
  if(sequential) ccd.info <- ccd.seq(ddatax,nc,m,alpha)
  # else ccd.info <- ccd.nonseq(ddatax,nc,m)
  else ccd.info <- ccd.nonseq_scale(ddatax,nc,m)
  
  r <- ccd.info$R   # radii of each object
  # M <- matrix(as.integer(ddatax < r+mac.eps), length(r))  ### label the covered objects with 0 & 1
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

# incrementally add cluster points to find the clustering with maximum silhouette
# old name ccd3.silhouette
# graph is the digraph CCD
# ddatax is the distance matrix
# cls is the actual classes 
# ind is the set of indices of clusterings to be checked
# lenDlimit is the maximum index of the dominating point to be checked for silhouette
# min.cls: the minimum percentage accepted as a cluster
ksccd.silhouette <- function(graph, ddatax, cls=NULL, min.cls=0, ind=NULL, lenDlimit=Inf){
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


# find the stongly connected components with given the radius of each object 
# old name ccd2.clustering
# datax, the data set
# m, K-S statistics parameter
# sequential, logical whether to use sequential K-S or not
ksccd.connected <- function(datax,m,sequential=FALSE,alpha=0.05){
  if(is.null(dim(datax))){
    return(list(member=1,R=0))
  } else {
  # info
  nr <- nrow(datax)
  nc <- ncol(datax)
  ddatax <- as.matrix(dist(datax))
  
  # find all radii 
  if(sequential) ccd.info <- ccd.seq(ddatax,nc,m,alpha)
  else ccd.info <- ccd.nonseq_scale(ddatax,nc,m)
  
  r <- ccd.info$R
  #M <- matrix(as.integer(ddatax < r+mac.eps), length(r)) #Matrix with entries 0,1 recording the object covered by each covering ball#
  M <- matrix(as.integer(ddatax <= r), length(r)) #Matrix with entries 0,1 recording the object covered by each covering ball#
  score <- rowSums(M)
  
  # find the intersection graph
  M.inter <- diag(T,nrow(M))
  for(i in 1:(nrow(M)-1)){
    for(j in (i+1):nrow(M)){
      temp <- (M[i,j] & M[j,i])
      if(temp) M.inter[i,j] <- M.inter[j,i] <- T
    }
  }
  
  # find the disconnected components that are the clusters
    member <- components.mat(M.inter)

  #R = r*dist.m
  return(list(member=member,R=r))
  }
}





# ccd clustering that find the connected components with greedy algorithm 
# old name ccd2.clustering
# datax, the data set
# m, K-S statistics parameter
# sequential, logical whether to use sequential K-S or not
# dom.method: 
#     greedy, dominating sets greedy on outdegree
#     ks    , dominating sets greedy on ks statistics
ksccd.clustering.connected <- function(datax,m,sequential=FALSE,
                                       dom.method="greedy2",alpha=0.05){
  # info
  nr <- nrow(datax)
  nc <- ncol(datax)
  ddatax <- as.matrix(dist(datax))
  dist.m <- max(ddatax)
  ddatax <- ddatax/dist.m
  mac.eps <- .Machine$double.eps
  
  # find all radii, and the dominating set 
  if(sequential) ccd.info <- ccd.seq(ddatax,nc,m,alpha)
  else ccd.info <- ccd.nonseq(ddatax,nc,m)
  
  r <- ccd.info$R
  M <- matrix(as.integer(ddatax < r+mac.eps), length(r)) #Matrix with entries 0,1 recording the object covered by each covering ball#
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
  R <- r[D] #D: the order of MDS#
  
  # find the intersection graph and its dominating set (M.dom)
  if(length(D)==1){
    D.ind <- 1
    D.member <- 1  # modified from NULL to 1 ##
  }
  else{
    M <- M[D,]
    M.dom <- diag(T,nrow(M))
    for(i in 1:(nrow(M.dom)-1)){
      for(j in (i+1):nrow(M.dom)){
        temp <- (M[i,] & M[j,])
        if(any(temp)) M.dom[i,j] <- M.dom[j,i] <- T
      }
    }
    # find the disconnected components that are the clusters
    D.member <- components.mat(M.dom)
  }
  R = R*dist.m
  return(list(D.member=D.member,D=D,R=R))
}


# ccd function that finds the radii of all points, sequential version
ccd.seq <- function(dx,d,m,alpha){
  
  n <- nrow(dx)
  p.scores <- 1:n
  ks <- sqrt(qchisq(1-alpha,2)/(4*n))
  
  res <- apply(dx,1,function(t){
    o <- order(t)
    r <- t[o]
    rw <- (p.scores/n)-m*(r/max(r))^d
    a <- match(TRUE,c(0,rw)[1:n]< -ks,nomatch=n)
    maxball <- which.max(rw[1:a])
    return(c(r[maxball],rw[maxball]))
  })
  res <- t(res)
  return(list(R=res[,1],KS=res[,2]))
}

# ccd function that finds the radii of all points, nonsequential version
# dx: distance matrix
# d: dimensionality
# m: KS null hypothesis parameter
ccd.nonseq <- function(dx,d,m){
  
  n <- nrow(dx)
  p.scores <- 1:n
  
  res <- apply(dx,1,function(t){
    o <- order(t)
    r <- t[o]
    rw <- (p.scores/n)-m*r^d
    maxball <- which.max(rw)
    return(c(r[maxball],rw[maxball]))
  })
  res <- t(res)
  return(list(R=res[,1],KS=res[,2]))
}


# ccd function that finds the radii of all points, nonsequential version
# dx: distance matrix
# d: dimensionality
# m: KS null hypothesis parameter
# no scaling!!!!!!!!!
ccd.nonseq_scale <- function(dx,d,m){
  
  n <- nrow(dx)
  p.scores <- 1:n
  
  res <- apply(dx,1,function(t){
    o <- order(t)
    r <- t[o]
    # rw <- p.scores-m*r^d
    rw <- p.scores-(m*r)^d
    maxball <- which.max(rw)
    return(c(r[maxball],rw[maxball]))
  })
  res <- t(res)
  return(list(R=res[,1],KS=res[,2]))
}



# which point is in which ball of dominators
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