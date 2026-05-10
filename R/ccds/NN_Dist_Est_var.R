library(MASS)

# data sets are simulated in a unit ball/sphere
# simulate a homogenuous poisson process in unit sphere
# requires library MASS
# d: dimension
# n: the size of data set
# critial quantile
rpoisball.unit <- function(n,d){
  # inner ball and outer ball values
  r1 <- runif(n,0,1)^(1/d)
  norm.data <- matrix(mvrnorm(n,rep(0,d),diag(d)),ncol=d,byrow=T)
  data1 <- apply(norm.data,1,function(x) x/sqrt(sum(x^2)))
  data1 <- apply(data1,1,function(x) x*r1)
  return(data1)
}

# Simulate the lower tail quantile of NN distance
NNDest.simpois.lower.quant = function(n, d, quant=0.99, niter=100, shape="sphere"){
  NN.dist.ave.mat = NULL
  NN.dist.var.mat = NULL
  
  if(shape=="sphere"){
    for(i in 1:niter){
      data.simu.list = lapply(1:n, rpoisball.unit, d=d)
      NN.dist.temp = sapply(2:n, function(x){
        data.temp = data.simu.list[[x]]
        data.dist = as.matrix(dist(data.temp))
        diag(data.dist) = Inf
        NN.dist.ttemp = apply(data.dist, 1, min) # the Nearest Neighbor distance for each points
        NN.dist.ttemp.ave = mean(NN.dist.ttemp)
        NN.dist.ttemp.var = var(NN.dist.ttemp)
        return(c(NN.dist.ttemp.ave, NN.dist.ttemp.var))
      })
      NN.dist.temp.ave = c(0, NN.dist.temp[1,])
      NN.dist.temp.var = c(0, NN.dist.temp[2,])
      NN.dist.ave.mat = rbind(NN.dist.ave.mat, NN.dist.temp.ave)
      NN.dist.var.mat = rbind(NN.dist.var.mat, NN.dist.temp.var)
    }
  }
  
  quant.ave.lower = sapply(1:n, function(x){
    quant = quantile(NN.dist.ave.mat[,x],1-quant)
  })
  names(quant.ave.lower)=NULL
  
  quant.var.upper = sapply(1:n, function(x){
    quant = quantile(NN.dist.var.mat[,x],quant)
  })
  names(quant.var.upper)=NULL
  
  return(list(average = quant.ave.lower, variance = quant.var.upper))
}

# Simulate the lower tail quantile of NN distance
# parallel computing
NNDestP.simpois.lower.quant = function(n, d, quant=0.99, niter=100, shape="sphere"){
  NN.dist.ave.mat = NULL
  NN.dist.var.mat = NULL
  
  if(shape=="sphere"){
    SimuOnce = function(){
      # simulate CSR within a unit sphere
      rpoisball.unit <- function(n,d){
        # inner ball and outer ball values
        r1 <- runif(n,0,1)^(1/d)
        norm.data <- matrix(mvrnorm(n,rep(0,d),diag(d)),ncol=d,byrow=T)
        data1 <- apply(norm.data,1,function(x) x/sqrt(sum(x^2)))
        data1 <- apply(data1,1,function(x) x*r1)
        return(data1)
      }
      
      data.simu.list = lapply(1:n, rpoisball.unit, d=d)
      NN.dist.temp = sapply(2:n, function(x){
        data.temp = data.simu.list[[x]]
        data.dist = as.matrix(dist(data.temp))
        diag(data.dist) = Inf
        NN.dist.ttemp = apply(data.dist, 1, min) # the Nearest Neighbor distance for each points
        NN.dist.ttemp.ave = mean(NN.dist.ttemp)
        NN.dist.ttemp.var = var(NN.dist.ttemp)
        return(c(NN.dist.ttemp.ave, NN.dist.ttemp.var))
      })
      NN.dist.temp.ave = c(0, NN.dist.temp[1,])
      NN.dist.temp.var = c(0, NN.dist.temp[2,])
      return(c(NN.dist.temp.ave=list(NN.dist.temp.ave), NN.dist.temp.var=list(NN.dist.temp.var)))
    }
    
    cores = detectCores()
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    NNest.m = foreach(1:niter,.packages = c("MASS","cluster","igraph")) %dopar% SimuOnce()
    stopCluster(cl)
    
    for(i in 1:niter){
      NN.dist.ave.mat = rbind(NN.dist.ave.mat, NNest.m[[i]]$NN.dist.temp.ave)
      NN.dist.var.mat = rbind(NN.dist.var.mat, NNest.m[[i]]$NN.dist.temp.var)
    }
  }
  
  quant.ave.lower = sapply(1:n, function(x){
    quant = quantile(NN.dist.ave.mat[,x],1-quant)
  })
  names(quant.ave.lower)=NULL
  
  quant.var.upper = sapply(1:n, function(x){
    quant = quantile(NN.dist.var.mat[,x], quant)
  })
  names(quant.var.upper)=NULL
  
  return(list(average = quant.ave.lower, variance = quant.var.upper))
}

# Simulate the upper tail quantile of NN distance
NNDest.simpois.upper.quant = function(n, d, quant=0.99, niter=100, shape="sphere"){
  NN.dist.ave.mat = NULL
  NN.dist.var.mat = NULL
  if(shape=="sphere"){
    for(i in 1:niter){
      data.simu.list = lapply(X=1:n, rpoisball.unit, d=d)
      set.seed(i)
      NN.dist.temp = sapply(2:n, function(x){
        data.temp = data.simu.list[[x]]
        data.dist = as.matrix(dist(data.temp))
        diag(data.dist) = Inf
        NN.dist.ttemp = apply(data.dist, 1, min) # the Nearest Neighbor distance for each points
        NN.dist.ttemp.ave = mean(NN.dist.ttemp)
        NN.dist.ttemp.var = variance(NN.dist.ttemp)
        return(c(NN.dist.ttemp.ave, NN.dist.ttemp.var))
      })
      NN.dist.temp.ave = c(0, NN.dist.temp[1,])
      NN.dist.temp.var = c(0, NN.dist.temp[2,])
      NN.dist.ave.mat = rbind(NN.dist.ave.mat, NN.dist.temp.ave)
      NN.dist.var.mat = rbind(NN.dist.var.mat, NN.dist.temp.var)
    }
  }
  
  quant.ave.upper = sapply(1:n, function(x){
    quant = quantile(NN.dist.ave.mat[,x], quant)
  })
  names(quant.ave.upper)=NULL
  
  quant.var.upper = sapply(1:n, function(x){
    quant = quantile(NN.dist.var.mat[,x],quant)
  })
  names(quant.var.upper)=NULL
  
  return(list(average = quant.ave.upper, variance = quant.var.upper))
}




# Simulate the two sided quantile of NN distance
NNDest.simpois.two.quant = function(n, d, quant=0.99, niter=100, shape="sphere"){
  NN.dist.ave.mat = NULL
  NN.dist.var.mat = NULL
  if(shape=="sphere"){
    for(i in 1:niter){
      data.simu.list = lapply(X=1:n, rpoisball.unit, d=d)
      set.seed(i)
      NN.dist.temp = sapply(2:n, function(x){
        data.temp = data.simu.list[[x]]
        data.dist = as.matrix(dist(data.temp))
        diag(data.dist) = Inf
        NN.dist.ttemp = apply(data.dist, 1, min) # the Nearest Neighbor distance for each points
        NN.dist.ttemp.ave = mean(NN.dist.ttemp)
        NN.dist.ttemp.var = var(NN.dist.ttemp)
        return(c(NN.dist.ttemp.ave, NN.dist.ttemp.var))
      })
      NN.dist.temp.ave = c(0, NN.dist.temp[1,])
      NN.dist.temp.var = c(0, NN.dist.temp[2,])
      NN.dist.ave.mat = rbind(NN.dist.ave.mat, NN.dist.temp.ave)
      NN.dist.var.mat = rbind(NN.dist.var.mat, NN.dist.temp.var)
    }
  }
  
  quant.ave.lower = sapply(1:n, function(x){
    quant = quantile(NN.dist.ave.mat[,x], (1-quant)/2)
  })
  names(quant.ave.lower)=NULL
  
  quant.ave.upper = sapply(1:n, function(x){
    quant = quantile(NN.dist.ave.mat[,x], (1+quant)/2)
  })
  names(quant.ave.upper)=NULL
  
  quant.var.lower = sapply(1:n, function(x){
    quant = quantile(NN.dist.var.mat[,x], (1-quant)/2)
  })
  names(quant.var.lower)=NULL
  
  quant.var.upper = sapply(1:n, function(x){
    quant = quantile(NN.dist.var.mat[,x], (1+quant)/2)
  })
  names(quant.var.upper)=NULL
  
  return(list(average.lower = quant.ave.lower,
              average.upper = quant.ave.upper,
              variance.lower = quant.var.lower,
              variance.upper = quant.var.upper))
}


# given a dist matrix, calculate the NN dist
# sc: scaling factor
# ddx: the dist matrix
NNDest.dist.f = function(ddx, sc){
  diag(ddx) = Inf
  dist.min = apply(ddx, 1, min) # the Nearest Neighbor distance for each points
  dist.min.ave =  mean(dist.min)/sc
  dist.min.var =  var(dist.min)/(sc^2)
  return(list(averge=dist.min.ave, variance=dist.min.var))
}

# Simulate the samples of NN distance
# n: maximum size of a data set
# d: dimension
# niter: number of data size to simulate
# shape: the shape of the unit support
NNDest.simpois = function(n, d, niter=100, shape="sphere"){
  if(shape=="sphere"){
    NN.dist.simu = lapply(X=1:n, niter=niter, function(X,niter){
      NN.dist.temp = NULL
      for(i in 1:niter){
        if(X==1){NN.dist.ttemp = 0; NN.dist.temp = c(NN.dist.temp, NN.dist.ttemp)}
        else {
          data.temp = rpoisball.unit(X, d)
          dist.mat = as.matrix(dist(data.temp))
          diag(dist.mat)=Inf
          NN.dist.ttemp = apply(dist.mat,1,min)
          NN.dist.temp = c(NN.dist.temp, NN.dist.ttemp)
        }
      }
      NN.dist = NN.dist.temp
      names(NN.dist) = NULL
      return(NN.dist)
    })
    return(NN.dist.simu)
  }
}
#NNDest.simpois(n=100,d=3,niter=10)

# Simulate the distance samples of a CSR
# n: maximum size of a data set
# d: dimension
# niter: number of data size to simulate
# shape: the shape of the unit support
dist.simpois = function(n, d, niter=10, shape="sphere"){
  if(shape=="sphere"){
    dist.simu = lapply(X=1:n, niter=niter, function(X,niter){
      dist.temp = NULL
      for(i in 1:niter){
        if(X==1){dist.ttemp = 0; dist.temp = c(dist.temp, dist.ttemp)}
        else {
          data.temp = rpoisball.unit(X, d)
          dist.mat = as.matrix(dist(data.temp))
          dist.ttemp = dist.mat[which(upper.tri(dist.mat))]
          dist.temp = c(dist.temp, dist.ttemp)
        }
      }
      dist = dist.temp
      return(dist)
    })
    return(dist.simu) 
  }
}
