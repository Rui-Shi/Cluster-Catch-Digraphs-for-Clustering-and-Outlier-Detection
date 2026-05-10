# auxiliary function for CCD clustering
library(igraph)
dominate <- function (data,distx,g,method = "greedy") 
{
  METHODS = c("greedy", "ks")
  method <- pmatch(tolower(method), METHODS)
  if (is.na(method)) {
    stop("invalid method")
  }
  if (method == 1) 
    dom <- dominate.greedy(data,g)
  return(x=dom$x,r=dom$r)
}

dominate.greedy <- function(g,weight=NULL)
{
  A <- g$A
  od <- apply(A,1,sum)+1
  S <- NULL
  diag(A) <- 0
  n <- nrow(A)
  covered <- rep(0,n)
  while(sum(covered)<n){
    i <- which.max(od)
    if(!is.null(weight)){
      qq <- weight
      qq[od != od[i]] <- -Inf
      i <- which.max(qq)
    }
    covered[A[i,]>0] <- 1
    covered[i] <- 1
    S <- c(S,i)
    A[,covered>0] <- 0
    od <- apply(A,1,sum)+1-covered
  }
  list(x=data[S,],r=g$R[S,])
}

# given the adjacency of digraph, 
# find the approximate minimum dominating set
# by means of the greedy algorithms
dominate.mat.greedy <- function(A)
{
  if(!is.matrix(A)) A <- matrix(A,nrow=1)
  S <- NULL
  n <- nrow(A)
  covered <- rep(FALSE,n)
  while(!all(covered)){
    od <- apply(A,1,sum)
    i <- which.max(od)
    covered[A[i,]==TRUE] <- TRUE
    S <- c(S,i)
    A[,covered==TRUE] <- FALSE
  }
  return(S)
}

# given the adjacency of digraph, 
# find the approximate minimum dominating set
# deleted points could be the set of dominating points
# take the maximum degree of deleted or non-deleted point
# by means of the greedy algorithms
dominate.mat.greedy2 <- function(A)
{
  if(!is.matrix(A)) A <- matrix(A,nrow=1)
  S <- NULL
  n <- nrow(A)
  covered <- rep(FALSE,n)
  score <- apply(A,1,sum)
  ind <- 1:n
  while(!all(covered)){
    if(is.null(S)){ 
      i <- which.max(score)
    } else {
      i <- ind2[which.max(score[-S])]
    }
    #covered[A[i,]==TRUE] <- TRUE
    covered[i] <- TRUE
    S <- c(S,i)
    ind2 <- ind[-S]
  }
  return(S)
}

# given the adjacency of digraph, 
# find the approximate minimum dominating set
# by means of the greedy algorithms
# deleted points could be the set of dominating points
# take only those that cover non-deleted 
dominate.mat.greedy3 <- function(A)
{
  if(!is.matrix(A)) A <- matrix(A,nrow=1)
  S <- NULL
  n <- nrow(A)
  covered <- rep(FALSE,n)
  while(!all(covered)){
    od <- apply(A,1,sum)
    i <- which.max(od)
    covered[A[i,]==TRUE] <- TRUE
    S <- c(S,i)
    A[,covered==TRUE] <- FALSE
  }
  return(S)
}

# find the approximate minimum dominating set 
# of a digraph with its adjacency matrix and its K-S info
dominate.mat.ks <- function(A,ks)
{
  if(!is.matrix(A)) A <- matrix(A,nrow=1)
  S <- NULL
  n <- nrow(A)
  covered <- rep(FALSE,n)
  while(!all(covered)){
    i <- which.max(ks)
    covered[A[i,]==TRUE] <- TRUE
    S <- c(S,i)
    A[,covered==TRUE] <- FALSE
    ks[covered==TRUE] <- -Inf
  }
  return(S)
}

# find disconnected sets of a graph given by the 
# adjacency matrix 
components.mat <- function(A){
  
  if(!is.matrix(A)) A <- matrix(A,nrow=1)
  
  # find the number of maximal connected parts of the graph
  g <- igraph::graph_from_adjacency_matrix(A)
  result <- igraph::components(g)$membership
  
  return(result)
}