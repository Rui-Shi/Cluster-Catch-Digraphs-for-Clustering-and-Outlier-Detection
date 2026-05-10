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
                            dom.method="greedy",alpha=0.05){
  
  # info
  nr <- nrow(datax)
  nc <- ncol(datax)
  ddatax <- as.matrix(dist(datax))
  ddatax <- ddatax/max(ddatax)
  mac.eps <- .Machine$double.eps   #########?????
  
  # find all radii, and the dominating set 
  if(sequential) ccd.info <- ccd.seq(ddatax,nc,m,alpha)
  else ccd.info <- ccd.nonseq(ddatax,nc,m)
  
  r <- ccd.info$R
  M <- matrix(as.integer(ddatax < r+mac.eps), length(r))
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
  if(length(D)==1) D.ind <- 1
  else{
    M <- M[D,]
    M.dom <- diag(T,nrow(M))
    for(i in 1:(nrow(M.dom)-1)){
      for(j in (i+1):nrow(M.dom)){
        M.dom[i,j] <- M.dom[j,i] <- any(M[i,] & M[j,])
      }
    }
    D.ind <- dominate.mat.ks(M.dom,score[D]) 
  }

  # keep the 2nd layer of dominating sets 
  Int.D=D[D.ind]
  Int.R=R[D.ind]
  
  # record scores 
  M <- matrix(M,nrow=length(D))
  score = rowSums(M)
  
  return(list(Int.D=Int.D,Int.R=Int.R,D=D,R=R,score=score))
}

# silhouette score based detection of the number of minimum dominating set
# incrementally add cluster points to find the clustering with maximum silhouette
# old name = ccd1.silhouette
ksccd.silhouette <- function(graph,ddatax,cls,ind=NULL,lenDlimit=NULL){
  
  lenD <- length(graph$score[graph$score>1])
  if(lenDlimit < lenD) lenD <- lenDlimit
  # print(graph$score[graph$score>1])
  dgraph <- graph
  if(is.null(ind)) ind <- 1:nrow(ddatax)
  
  maxsi <- 0
  si.ind <- 0
  
  for(i in 2:lenD){
    dgraph$Int.D <- graph$Int.D[1:i]
    dgraph$Int.R <- graph$Int.R[1:i]
    result <- ccd3.clustering.nonvalid(dgraph,ddatax)
    if(length(unique(result[ind])) < 2) datasi <- 0
    else{
      datasi <- silhouette(result[ind],ddatax[ind,ind])
      datasi <- mean(datasi[,3])  
    }
    if(datasi > maxsi){
      maxsi <- datasi
      si.ind <- i
    } 
    print(paste("sil", i," = ", datasi,sep=""))
  }
  
  return(list(si=maxsi,si.ind=si.ind))
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
  mac.eps <- .Machine$double.eps
  
  # find all radii, and the dominating set 
  if(sequential) ccd.info <- ccd.seq(ddatax,nc,m,alpha)
  else ccd.info <- ccd.nonseq(ddatax,nc,m)
  
  r <- ccd.info$R
  M <- matrix(as.integer(ddatax < r+mac.eps), length(r))
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
  if(length(D)==1){
    D.ind <- 1
    D.member <- NULL
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
  
  
  return(list(D.member=D.member,D=D,R=R))
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