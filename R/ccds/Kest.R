library(MASS)
# auxiliary functions for simulations and Ripley's K function

# simulate a homogenuous poisson process in unit sphere
# requires library MASS
rpoisball.unit <- function(n,d){
  # inner ball and outer ball values
  r1 <- runif(n,0,1)^(1/d)
  norm.data <- matrix(mvrnorm(n,rep(0,d),diag(d)),ncol=d,byrow=T)
  data1 <- apply(norm.data,1,function(x) x/sqrt(sum(x^2)))
  data1 <- apply(data1,1,function(x) x*r1)
  return(data1)
}

# simulate a homogenuous poisson process in unit sphere
# normal distribution method, requires MASS library
# this is a version, do not need any edge correction
# requires library MASS
# m is number of obs, d is number of dimensions
rpoisball <- function(n,d){
  
  # check properties for dimensionality d
  m <- (2^d-1)*n
  
  # inner ball and outer ball values
  r1 <- runif(n,0,1)^(1/d)
  r2 <- (runif(m,0,1)^(1/d)+1)
  
  norm.data <- mvrnorm(n+m,rep(0,d),diag(d))
  data1 <- apply(norm.data[1:n,],1,function(x) x/sqrt(sum(x^2)))
  data2 <- apply(norm.data[(n+1):(n+m),],1,function(x) x/sqrt(sum(x^2)))
  data1 <- apply(data1,1,function(x) x*r1)
  data2 <- apply(data2,1,function(x) x*r2)
  
  return(rbind(data1,data2))
}

# simulate a homogenuous poisson process in unit box
# requires library MASS
rpoisbox.unit <- function(n,d){
  
  # independent uniform distributions for each dimension
  data1 <- sapply(1:d, function(x) runif(n,0,1))

  return(data1)
}

# simulate points in a box with at least some distance apart
# requires library: flexclust
# n: number of points
# r: at least distance between points
# d: dimensions
rUnifDist <- function(n,r,d,method){
  
  datax <- NULL
  count <- 0
  repeat{
    count <- count + 1
    
    # if null, generate a point, else check if the point is at least r distant
    if(is.null(datax)){
      datax <- matrix(runif(d),nrow=1)
    } else {
      datay <- matrix(runif(d),nrow=1)
      dxy <- dist2(datax,datay,method)
      if(all(dxy > r)) datax <- rbind(datax,datay)
    }
    if(nrow(datax)> (n-1)) break

    # if cant find any point r distant in 100 trials, reset
    if(count==100){
      datax <- NULL
      count <- 0
    }
  }
  
  return(datax)
}

# simulate Box Matern process in a unit box
# all points are in the window, it is a special kind of matern process
# nclus: number of clusters
# mu: # of obs in cluster
# r: radius of cluster
# d: dimension
rMaternUnif <- function(nclus,mu,r,rperc,d,dist.method="maximum"){
  
  # random uniform points in unit box, with at least 0.5*r apart
  dataclus <- rUnifDist(nclus,(2+rperc)*r,d,method=dist.method)

  # blocks of random boxes, scale to r
  # databall <- rpoisball.unit(mu*nclus,d)*r
  databall <- matrix(runif(mu*nclus*d,-r,r),ncol=d)
  
  # scale to r, adjust to centers and record
  datax <- NULL
  for(i in 1:nclus){
    temp <- databall[(1:mu)+(i-1)*mu,]
    temp <- t(apply(temp,1,function(x) x+dataclus[i,]))
    datax <- rbind(datax,temp)  
  }
  
  return(datax)
}

# simulate Box Matern process in a unit box, each is cluster is normal
# all points are in the window, it is a special kind of matern process
# nclus: number of clusters
# mu: # of obs in cluster
# r: radius of cluster
# d: dimension
rMaternNormal <- function(nclus,mu,r,rperc,d,dist.method="maximum"){
  
  # random uniform points in unit box, with at least 0.5*r apart
  dataclus <- rUnifDist(nclus,(2+rperc)*r,d,method=dist.method)
  
  # blocks of random boxes, scale to r
  # databall <- rpoisball.unit(mu*nclus,d)*r
  databall <- mvrnorm(mu*nclus,rep(0,d),diag(r/20,d))
  
  # scale to r, adjust to centers and record
  datax <- NULL
  for(i in 1:nclus){
    temp <- databall[(1:mu)+(i-1)*mu,]
    temp <- t(apply(temp,1,function(x) x+dataclus[i,]))
    datax <- rbind(datax,temp)  
  }
  
  return(datax)
}

# putting random noise in between clusters
# datax: unlabelled data
# cls: class information
# n: number of points to be simulated as noise
rNoise <- function(datax,cls,n){
  
  # for each class divide the number of noisy points
  uni.cls <- unique(cls)
  nc <- length(uni.cls)
  comb.n <- combn(nc,2)
  ncol.comb <- ncol(comb.n)
  n <- round(n/ncol.comb)
  
  # random pairs of noise 
  r.ind <- sample(1:(nrow(datax)/nc),n)
  
  # add random points between points of different clusters
  datan <- NULL
  for(i in 1:ncol.comb){
    t <- comb.n[,i]
    data1 <- datax[which(cls==uni.cls[t[1]]),]
    data2 <- datax[which(cls==uni.cls[t[2]]),]
    data1 <- data1[r.ind,]
    data2 <- data2[r.ind,]
    theta <- runif(n)
    temp <- data1*theta+data2*(1-theta)
    datan <- rbind(datan,temp)
  }

  return(datan)
}

# putting random noise between booundaries of clusters
# datax: unlabelled data
# cls: class information
# n: number of points to be simulated as noise
# sigma: sigma parameter of the noisy gauss dist
rNoise.bound <-   function(datax,cls,n){
  
  # for each class divide the number of noisy points
  uni.cls <- unique(cls)
  nc <- length(uni.cls)
  d <- ncol(datax)

  # get n number of boundary points from each clusters,
  # given the distances to center of the clusters
  data.bound <- matrix(nrow=nc,ncol=n)
  for(i in 1:nc){
    ind <- which(cls==uni.cls[i])
    datan <- datax[ind,]
    data.mean <- matrix(apply(datan,2,mean),nrow=1)
    dist.data <- dist2(datan,data.mean)
    data.bound[i,] <- ind[order(dist.data)[1:n]] 
  }
  
  # given the boundary points, generate random points
  # between these points
  datan <- apply(data.bound,2,function(x){
    p.bound <- datax[x,]
    p.rand <- apply(p.bound,2,function(t){
      r.coef <- runif(nc)
      r.coef <- r.coef/sum(r.coef)
      return(sum(t*r.coef))                                                    
    })
  })
  datan <- t(datan)
  
  return(datan)
  
}

# putting random noise on centers of clusters, Gaussian noise
# datax: unlabelled data
# cls: class information
# n: number of points to be simulated as noise
# sigma: sigma parameter of the noisy gauss dist
rNoise.addGauss <-   function(datax,cls,n,sigma){
  
  # for each class divide the number of noisy points
  uni.cls <- unique(cls)
  nc <- length(uni.cls)
  d <- ncol(datax)
  
  # for each clusters center, add a noise
  datanew <- NULL
  for(i in 1:nc){
    temp <- datax[which(cls==uni.cls[i]),]
    meantemp <- apply(temp,2,mean)
    datan <- mvrnorm(n,meantemp,diag(sigma,d))
    datanew <- rbind(datanew,datan) 
  }
  
  return(datanew)
}

# add Noise to the entire domain
# n: number of points
# r: at least distance between points
# d: dimensions
rNoise.entire <- function(n,a,b,d){
  
  datax <- matrix(runif(n*d,a,b),ncol=d)
  
  return(datax)
}

# simulate a strauss process in unit box 
strauss<-function(n,r,g,d){
  repeat {
    kappa <- n
    mypar <- list(beta = kappa, r=r, gamma = g)
    mo <- list(cif = "strauss", par = mypar, w = window)
    mypp <- rmh(model = mo, start = list(n.start = n), control = list(nrep = 100000, p=1))
    if(mypp$n>0){
      write.table(mypp,file="mypp.txt",col.names=T,row.names=FALSE)
      temp<-read.table("mypp.txt",header=T)
    }else{
      next
    }
    if (nrow(temp) ==n) break
  }
  return(mypp)
  
} #end of strauss

# get the Kest envelopes provided by the simulations, new version 2
# no edge correction, correction commented out
# (the slopes of envelopes)
# m is the number of observations in the original data set
# rn is the number of lengths in radii of Kest
# ln is the smallest cardinality ball (dont use it for now)
# niter is the number of iterations
Kest.simpois <- function(m,d,rn,niter){
  
  # do this for each simulated data set
  # record points for simulated cases
  r <- seq(1/rn,1,1/rn)
  c <- (2^d-1)
  mat.flag <- sapply(1:m,function(t) c(rep(FALSE,c*t),rep(TRUE,c*(m-t))),simplify=TRUE)
  mat.flag2 <- sapply(1:(m-1),function(t) c(rep(FALSE,c*t),rep(TRUE,c),rep(FALSE,c*(m-t-1))),simplify=TRUE)
  
  # simulation
  Kest.m <- NULL
  for(i in 1:niter){
    
    # simulating data
    # print(i)
    temp <- rpoisball.unit(m,d)
    
    # distances
    temp.dist1 <- as.matrix(dist(temp))
    diag(temp.dist1) <- Inf

    # analyze
    result <- sapply(r,function(x){
      Mtemp1 <- (temp.dist1 < x)
      Mtemp1[lower.tri(Mtemp1)] <- 0
      sumM1 <- cumsum(2*colSums(Mtemp1))
      return(sumM1/((1:m)*(1:m)))
    },simplify=TRUE) # simple = T: return a matrix rather than list.
    Kest.m <- rbind(Kest.m,as.vector(result))
  }
  
  Kest.max <- apply(Kest.m,2,max)
  Kest.min <- apply(Kest.m,2,min)
  Kest.max <- matrix(Kest.max,nrow=m)
  Kest.min <- matrix(Kest.min,nrow=m)
  
  return(list(max=Kest.max,min=Kest.min, r=r))
}

# get the Kest envelopes provided by the simulations, new version 2
# edge correction, translation correction
# (the slopes of envelopes)
# old name = Kest.simpois2
# m is the number of observations in the original data set
# rn is the number of lengths in radii of Kest
# ln is the smallest cardinality ball (dont use it for now)
# niter is the number of iterations
Kest.simpois.edge <- function(m,d,rn,niter){
  
  # do this for each simulated data set
  # record points for simulated cases
  r <- seq(1/rn,1,1/rn)

  # simulation
  Kest.m <- NULL  
  for(i in 1:niter){
    
    # simulating data
    # print(i)
    temp <- rpoisball.unit(m,d)
    
    # distances
    temp.dist <- as.matrix(dist(temp))

    # calculate weights for correction
    cons <- (sqrt(pi)*gamma((d+1)/2))/(2*gamma(d/2+1))
    integrand <- function(t) sin(t)^d
    ftemp <- sapply(temp.dist,function(t){
      return(integrate(integrand,0,acos(t/2))$value)
    },simplify = TRUE)
    ftemp <- matrix(ftemp,nrow=nrow(temp.dist),byrow = FALSE)
    ftemp <- cons*(1/ftemp)
    
    # analyze
    diag(temp.dist) <- Inf
    result <- sapply(r,function(x){
      Mtemp <- (temp.dist < x)
      Mtemp[lower.tri(Mtemp)] <- 0
      ftemp[lower.tri(ftemp)] <- 0
      Mtemp <- Mtemp*ftemp
      sumM <- cumsum(2*colSums(Mtemp))
      return(sumM/((1:m)*(1:m)))
    },simplify=TRUE)
    Kest.m <- rbind(Kest.m,as.vector(result))
  }
  
  Kest.max <- apply(Kest.m,2,max)
  Kest.min <- apply(Kest.m,2,min)
  Kest.max <- matrix(Kest.max,nrow=m)
  Kest.min <- matrix(Kest.min,nrow=m)
  
  return(list(max=Kest.max,min=Kest.min, r=r))
}

# get the Kest envelopes provided by the simulations, new version 2
# edge correction, translation correction
# (the slopes of envelopes)
# added support for quantiles given in simulated envelopes
# old name = Kest.simpois2.quantile
# m is the number of observations in the original data set
# rn is the number of lengths in radii of Kest
# quan is the set of quantiles you wanna use
# niter is the number of iterations
Kest.simpois.edge.quantile <- function(m,d,rn,quan,niter){

  # do this for each simulated data set
  # record points for simulated cases

  r <- seq(1/rn,1,1/rn)
  # simulation
  Kest.m <- NULL
  for(i in 1:niter){
    # i = 10
    # simulating data
    #print(i)
    temp <- rpoisball.unit(m,d)

    # distances
    temp.dist <- as.matrix(dist(temp))

    # calculate weights for correction
    cons <- (sqrt(pi)*gamma((d+1)/2))/(2*gamma(d/2+1))
    integrand <- function(t) sin(t)^d
    ftemp <- sapply(temp.dist,function(t){
      return(integrate(integrand,0,acos(t/2))$value)
    },simplify = TRUE)
    ftemp <- matrix(ftemp,nrow=nrow(temp.dist),byrow = FALSE)
    ftemp <- cons*(1/ftemp)

    # analyze
    diag(temp.dist) <- Inf
    result <- sapply(r,function(x){
      Mtemp <- (temp.dist < x)
      Mtemp[lower.tri(Mtemp)] <- 0
      ftemp[lower.tri(ftemp)] <- 0
      Mtemp <- Mtemp*ftemp
      sumM <- cumsum(2*colSums(Mtemp))
      return(sumM/((1:m)*(1:m)))
    },simplify=TRUE)
    Kest.m <- rbind(Kest.m,as.vector(result))
  }

  Kest.quan <- list()
  for(cur_quan in quan){
    temp <- apply(Kest.m,2, quantile, probs = as.numeric(cur_quan))
    Kest.quan[[as.character(cur_quan)]] <- matrix(temp,nrow=m)
  }

  return(list(Kest.m = Kest.m, quan=Kest.quan, r=r))
}

# get the Kest envelopes provided by the simulations, new version 2
# edge correction, translation correction
# (the slopes of envelopes)
# added support for quantiles given in simulated envelopes
# old name = Kest.simpois2.quantile 
# m is the number of observations in the original data set
# rn is the number of lengths in radii of Kest
# quan is the set of quantiles you wanna use
# niter is the number of iterations
# parallel computing 
KestP.simpois.edge.quantile <- function(m,d,rn,quan,niter){
  
  # do this for each simulated data set
  # record points for simulated cases
  
  r <- seq(1/rn,1,1/rn)
  # simulation
  Kest.m <- NULL
  SimuOnce = function(){
    
    rpoisball.unit <- function(n,d){
      # inner ball and outer ball values
      r1 <- runif(n,0,1)^(1/d)
      norm.data <- matrix(mvrnorm(n,rep(0,d),diag(d)),ncol=d,byrow=T)
      data1 <- apply(norm.data,1,function(x) x/sqrt(sum(x^2)))
      data1 <- apply(data1,1,function(x) x*r1)
      return(data1)
    }
    temp <- rpoisball.unit(m,d)
    
    # distances
    temp.dist <- as.matrix(dist(temp))
    
    # calculate weights for correction
    cons <- (sqrt(pi)*gamma((d+1)/2))/(2*gamma(d/2+1))
    integrand <- function(t) sin(t)^d
    ftemp <- sapply(temp.dist,function(t){
      return(integrate(integrand,0,acos(t/2))$value)
    },simplify = TRUE)
    ftemp <- matrix(ftemp,nrow=nrow(temp.dist),byrow = FALSE)
    ftemp <- cons*(1/ftemp)
    
    # analyze
    diag(temp.dist) <- Inf
    result <- sapply(r,function(x){
      Mtemp <- (temp.dist < x)
      Mtemp[lower.tri(Mtemp)] <- 0
      ftemp[lower.tri(ftemp)] <- 0
      Mtemp <- Mtemp*ftemp
      sumM <- cumsum(2*colSums(Mtemp))
      return(sumM/((1:m)*(1:m)))
    },simplify=TRUE)
    
    return(as.vector(result))
  }
  cores = detectCores()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  Kest.m = foreach(1:niter,.packages = c("MASS","cluster","igraph")) %dopar% SimuOnce()
  stopCluster(cl)
  
  Kest.m <- do.call(rbind, Kest.m)

  Kest.quan <- list()
  for(cur_quan in quan){
  temp <- apply(Kest.m,2, quantile, probs = as.numeric(cur_quan))
  Kest.quan[[as.character(cur_quan)]] <- matrix(temp,nrow=m)
  }

  return(list(Kest.m = Kest.m, quan=Kest.quan, r=r))
}

# get the Kest envelopes provided by the simulations, new version 3
# the t for the Kest envelopes are dynamic, adapted for higher dimensions
# edge correction, translation correction
# (the slopes of envelopes)
# added support for quantiles given in simulated envelopes
# m is the number of observations in the original data set
# quan is the set of quantiles you wanna use
# niter is the number of iterations
Kest.simpois.edge.quantile.dynamic <- function(m,d,rn,quant,niter){
  #m=100;d=20;rn=5;quant=0.99;niter=100
  # do this for each simulated data set
  # record points for simulated cases
  data.temp = lapply(1:niter, function(x){
    simu = rpoisball.unit(m,d)
    return(simu)})
  
  # simulate the sample of distance between two points
  dist.simu = sapply(1:niter, function(x){
    dist.mat = as.matrix(dist(data.temp[[x]]))
    dist.temp = dist.mat[which(upper.tri(dist.mat))]
    return(dist.temp)
  })
  dist.simu = c(dist.simu)
  r = quantile(dist.simu, seq(0.5/rn,0.5,0.5/rn))
  
  # simulation
  Kest.m <- NULL
  for(i in 1:niter){
    # simulating data
    #print(i)
    temp <- data.temp[[i]]
    
    # distances
    temp.dist <- as.matrix(dist(temp))
    
    # calculate weights for correction
    cons <- (sqrt(pi)*gamma((d+1)/2))/(2*gamma(d/2+1))
    integrand <- function(t) sin(t)^d
    ftemp <- sapply(temp.dist,function(t){
      return(integrate(integrand,0,acos(t/2))$value)
    },simplify = TRUE)
    ftemp <- matrix(ftemp,nrow=nrow(temp.dist),byrow = FALSE)
    ftemp <- cons*(1/ftemp)
    
    # analyze
    diag(temp.dist) <- Inf
    result <- sapply(r,function(x){
      Mtemp <- (temp.dist < x)
      Mtemp[lower.tri(Mtemp)] <- 0
      ftemp[lower.tri(ftemp)] <- 0
      Mtemp <- Mtemp*ftemp
      sumM <- cumsum(2*colSums(Mtemp))
      return(sumM/((1:m)*(1:m)))
    },simplify=TRUE)
    Kest.m <- rbind(Kest.m,as.vector(result))
  }
  
  Kest.quan <- list()
  for(cur_quan in quant){
    temp <- apply(Kest.m,2, quantile, probs = as.numeric(cur_quan))
    Kest.quan[[as.character(cur_quan)]] <- matrix(temp,nrow=m)
  }
  
  return(list(Kest.m = Kest.m, quant=Kest.quan, r=r))
}

#test = Kest.simpois.edge.quantile.dynamic(m = 100, d=20, quant=0.99, niter = 100, rn = 5)
#test1 = Kest.simpois.edge.quantile(m = 100, d=20, quan=0.99, niter = 100, rn = 10)


# merge Kest simpois quantiles from multiple runs or objects
# each Kest should be an element of a list
Kest.simpois.edge.quantile.merge <- function(list_of_Kest,m,d,rn,quan,niter){
  
  # do this for each simulated data set
  # record points for simulated cases
  r <- seq(1/rn,1,1/rn)
  
  # simulation
  Kest.m <- NUL
  for(i in 1:length(list_of_Kest)){
    
    # simulating data
    print(i)
    Kest.m <- rbind(Kest.m,list_of_Kest[[i]])
  }
  
  # get Kest
  Kest.quan <- list()
  for(cur_quan in quan){
    temp <- apply(Kest.m,2, quantile, probs = as.numeric(cur_quan))
    Kest.quan[[as.character(cur_quan)]] <- matrix(temp,nrow=m)
  }
  
  return(list(Kest.m = Kest.m, quan=Kest.quan, r=r))
}

# get the Kest envelopes provided by the simulations for unit ball in a box
# no edge correction
# the unit box is of of (0,1)^d, with center (0.5, 0.5) and ball radious of  
# envelope distances are between 0 to 0.5 with 0.5/rn increments
# (the slopes of envelopes)
# m is the number of observations in the original data set
# rn is the number of lengths in radii of Kest
# ln is the smallest cardinality ball (dont use it for now)
# niter is the number of iterations
Kest.simpois.ballinbox <- function(m,d,rn,niter){
  
  # do this for each simulated data set
  # record points for simulated cases
  r <- seq(1/rn,1,1/rn)
  
  # simulation
  Kest.m <- NULL
  for(i in 1:niter){
    
    # simulating data until m number of points fall within the 
    # unit ball in a unit box
    # print(i)
    temp <- matrix(NA, nrow = 1, ncol = d)
    j <- 0
    while(nrow(temp) < m+1){
      temp_cur <- rpoisbox.unit(m,d)
      temp_cur <- (temp_cur*2)-1 # adjust box to (-1,1)^d for a unit ball 
      ind_inball <- apply(temp_cur, 1, function(x) return(sqrt(sum(x^2))))
      temp_cur <- temp_cur[ind_inball < 1,]
      temp <- rbind(temp, temp_cur)
      j <- j + 1
      # print(c(j, nrow(temp)))
    }
    temp <- na.omit(temp)
    temp <- temp[1:m,]
    
    # distances
    temp.dist <- as.matrix(dist(temp))
    diag(temp.dist) <- Inf
    
    # analyze
    result <- sapply(r,function(x){
      Mtemp1 <- (temp.dist < x)
      Mtemp1[lower.tri(Mtemp1)] <- 0
      sumM1 <- cumsum(2*colSums(Mtemp1))
      return(sumM1/((1:m)*(1:m)))
    },simplify=TRUE)
    Kest.m <- rbind(Kest.m,as.vector(result))
  }
  
  Kest.max <- apply(Kest.m,2,max)
  Kest.min <- apply(Kest.m,2,min)
  Kest.max <- matrix(Kest.max,nrow=m)
  Kest.min <- matrix(Kest.min,nrow=m)
  
  return(list(max=Kest.max,min=Kest.min, r=r))
}

# get the Kest envelopes provided by the simulations for unit ball in a box
# edge correction, translation correction
# the unit box is of of (0,1)^d, with center (0.5, 0.5) and ball radious of  
# envelope distances are between 0 to 0.5 with 0.5/rn increments
# (the slopes of envelopes)
# m is the number of observations in the original data set
# rn is the number of lengths in radii of Kest
# ln is the smallest cardinality ball (dont use it for now)
# niter is the number of iterations
Kest.simpois.ballinbox.withcorrection <- function(m,d,rn,niter){
  
  # do this for each simulated data set
  # record points for simulated cases
  r <- seq(1/rn,1,1/rn)
  
  # simulation
  Kest.m <- NULL
  for(i in 1:niter){
      
    # simulating data until m number of points fall within the 
    # unit ball in a unit box
    # print(i)
    temp <- matrix(NA, nrow = 1, ncol = d)
    j <- 0
    while(nrow(temp) < m+1){
      temp_cur <- rpoisbox.unit(m,d)
      temp_cur <- (temp_cur*2)-1 # adjust box to (-1,1)^d for a unit ball  
      ind_inball <- apply(temp_cur, 1, function(x) return(sqrt(sum(x^2))))
      temp_cur <- temp_cur[ind_inball < 1,]
      temp <- rbind(temp, temp_cur)
      j <- j + 1
      # print(c(j, nrow(temp)))
    }
    temp <- na.omit(temp)
    temp <- temp[1:m,]
    
    # distances
    temp.dist <- as.matrix(dist(temp))
    
    # calculate weights for correction
    cons <- (sqrt(pi)*gamma((d+1)/2))/(2*gamma(d/2+1))
    integrand <- function(t) sin(t)^d
    ftemp <- sapply(temp.dist,function(t){
      return(integrate(integrand,0,acos(t/2))$value)
    },simplify = TRUE)
    ftemp <- matrix(ftemp,nrow=nrow(temp.dist),byrow = FALSE)
    ftemp <- cons*(1/ftemp)
    
    # analyze
    diag(temp.dist) <- Inf
    result <- sapply(r,function(x){
      Mtemp <- (temp.dist < x)
      Mtemp[lower.tri(Mtemp)] <- 0
      ftemp[lower.tri(ftemp)] <- 0
      Mtemp <- Mtemp*ftemp
      sumM <- cumsum(2*colSums(Mtemp))
      return(sumM/((1:m)*(1:m)))
    },simplify=TRUE)
    Kest.m <- rbind(Kest.m,as.vector(result))
  }
  
  Kest.max <- apply(Kest.m,2,max)
  Kest.min <- apply(Kest.m,2,min)
  Kest.max <- matrix(Kest.max,nrow=m)
  Kest.min <- matrix(Kest.min,nrow=m)
  
  return(list(max=Kest.max,min=Kest.min, r=r))
}

# get the Kest envelopes provided by the simulations for unit box
# no edge correction
# (the slopes of envelopes)
# m is the number of observations in the original data set
# rn is the number of lengths in radii of Kest
# ln is the smallest cardinality ball (dont use it for now)
# niter is the number of iterations
Kest.simpois.box <- function(m,d,rn,niter){
  
  # do this for each simulated data set
  # record points for simulated cases
  r <- seq(1/rn,1,1/rn)
  
  # simulation
  Kest.m <- NULL
  for(i in 1:niter){
    
    # simulating data
    # print(i)
    temp <- rpoisbox.unit(m,d)
    
    # distances
    temp.dist <- as.matrix(dist(temp))
    diag(temp.dist) <- Inf
    
    # analyze
    result <- sapply(r,function(x){
      Mtemp1 <- (temp.dist < x)
      Mtemp1[lower.tri(Mtemp1)] <- 0
      sumM1 <- cumsum(2*colSums(Mtemp1))
      return(sumM1/((1:m)*(1:m)))
    },simplify=TRUE)
    Kest.m <- rbind(Kest.m,as.vector(result))
  }
  
  Kest.max <- apply(Kest.m,2,max)
  Kest.min <- apply(Kest.m,2,min)
  Kest.max <- matrix(Kest.max,nrow=m)
  Kest.min <- matrix(Kest.min,nrow=m)
  
  return(list(max=Kest.max,min=Kest.min, r=r))
}

# get the Kest envelopes provided by the simulations, new version
# (the slopes of envelopes)
# m is the number of observations in the original data set
# c is the chunk size for the simulation
# rn is the number of lengths in radii of Kest
# ln is the smallest cardinality ball (dont use it for now)
Kest.simpois.old2 <- function(m,d,c,rn,niter){
  
  # do this for each simulated data set
  r <- seq(1/rn,1,1/rn)

  Kest.m <- NULL
  for(i in 1:niter){
    # print(i)
    temp <- rpoisball.old(m,d,c)
    temp.dist <- as.matrix(dist(temp))
    diag(temp.dist) <- Inf
    M.temp.dist <- sapply(temp.dist,function(x){
      sumr <- sum(x < r)
      return(sumr)
    },simplify=TRUE)
    M.temp.dist <- matrix(M.temp.dist,ncol=ncol(temp.dist))
    Row.temp.dist <- t(apply(M.temp.dist,1,tabulate,nbins=length(r)))
    # print(Row.temp.dist)
    stop()
    Col.temp.dist <- apply(Row.temp.dist,2,cumsum)
    Col.temp.dist <- Col.temp.dist#/((1:m)*(0:(m-1)))
    Kest.m <- rbind(Kest.m,as.vector(Col.temp.dist))
  }
  
  Kest.max <- apply(Kest.m,2,max)
  Kest.min <- apply(Kest.m,2,min)
  Kest.max <- matrix(Kest.max,nrow=m)
  Kest.min <- matrix(Kest.min,nrow=m)
  
  return(list(max=Kest.max,min=Kest.min))
}

# get the Kest envelopes provided by the simulations, older version
# (the slopes of envelopes)
# sim.pois is the simulation results
# m is the number of observations in the original data set
# rn is the number of lengths in radii of Kest
# ln is the smallest cardinality ball (dont use it for now)
Kest.simpois.old <- function(sim.pois,m,rn){
  
  data <- sim.pois$data
  dist.data <- sim.pois$dist.data
  d <- sim.pois$dim
  br <- ((1:m)^(1/d))/(m^(1/d))
  r <- sapply(br,function(x) seq(x/rn,x,x/rn),simplify = TRUE)
  
  # do this for each simulated data set
  Kest.m <- replicate(m,NULL)
  for(i in 1:length(data)){
    temp <- data[[i]]
    temp.dist <- dist.data[[i]]
    
    # get the distances from the center and sort
    dist.cent <- apply(temp,1,function(x) sqrt(sum(x^2)))
    
    # for each poisson dist, find the curves
    for(j in 1:m){
      ind <- which(dist.cent < br[j])
      
      # if window has less than 2 elements, abort. 
      if(length(ind)<2) result <- rep(0,rn)
      else result <- Kest.f(temp.dist[ind,],r[,j],d)
      
      Kest.m[[j]] <- rbind(Kest.m[[j]],result)
    }
  }
  
  # find high and low curves
  slopes <- list()
  for(i in 1:m) slopes[[i]] <- apply(Kest.m[[i]],2,range)
  
  return(slopes)
}

# the Kest function
# adjust neighbor counts and adjust to number of points
# no edge correction
# M is the logical matrix that rows are points inside window,
#      columns are points in the bigger circle
# r is the radii of the window 
Kest.f <- function(M,r){
  
  # count self catching points
  n <- nrow(M)
  nc <- ncol(M)
  result <- sapply(r,function(x){
    temp <- (M<x)
    diag(temp) <- FALSE
    temp <- sum(temp)
    return(temp)
  },simplify = TRUE)
  result <- result/(n*nc)

  return(result)
}

# the Kest function
# just neighbor counts and adjust to number of points
# edge corrected version
# old name = Kest.f2
# M is the logical matrix that rows are points inside window,
#      columns are points in the bigger circle
# r is the radii of the window
# d is the dimensionality of the data
# sc is the scaling factor
Kest.f.edge <- function(M,r,sc,d){
  
  # get the edge correction weights
  cons <- (sqrt(pi)*gamma((d+1)/2))/(2*gamma(d/2+1))
  rmax <- max(r)
  #Mf <- M/rmax
  Mf <- M/sc

  integrand <- function(t) sin(t)^d
  Mf <- sapply(Mf,function(t){
    return(integrate(integrand,0,acos(t/2))$value)
  },simplify = TRUE)
  Mf <- matrix(Mf,nrow=nrow(M),byrow = FALSE)
  Mf <- cons*(1/Mf)
  Mf[Mf==Inf] <- max(Mf[Mf!=Inf])

  # count self catching points
  n <- nrow(M)
  nc <- ncol(M)
  result <- sapply(r,function(x){
    temp <- (M<x)
    diag(temp) <- FALSE
    temp <- sum(temp*Mf)
    return(temp)
  },simplify = TRUE)
  result <- result/(n*nc)
  
  
  # if you want to add area and number info,
  # add some of the codes below, but the codes above are enough for 
  # comparison
  #r <- r/max(r)
  #result <- ((pi*r^2)/(n^2-n))*result-pi*r^2

  return(result)
}