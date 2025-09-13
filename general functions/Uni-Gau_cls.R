library(spatstat)
library(hypercube)
library(MASS)
#library(special)

# simulate a homogenuous poisson process in unit sphere
# requires library MASS
rpoisball_unit <- function(n,d){
  # inner ball and outer ball values
  r1 <- runif(n,0,1)^(1/d)
  norm.data <- matrix(mvrnorm(n,rep(0,d),diag(d)),ncol=d,byrow=T)
  data1 <- apply(norm.data,1,function(x) x/sqrt(sum(x^2)))
  data1 <- apply(data1,1,function(x) x*r1)
  return(data1)
}

# the volume of a n-dim sphere
volume_n_sphere <- function(r, d) {
  pi_power <- pi^(d/2)
  gamma_value <- gamma(d/2 + 1)
  r_power <- r^d
  V <- (pi_power / gamma_value) * r_power
  return(V)
}

# simulate a data set with uniform clusters and noise (no outliers)
# n_cls: the number of each cluster
# d: dimension
# center: the center of each cluster by row
# rad_range: the upper & lower bound of random radius
# noise: whether there is noise
# r_noise: the radius of noise cluster
# n_noise: the number of noise
uniform_cls = function(n_cls=c(100,100), d=2, 
                       centers=rbind(c(0, 0), c(0, 3)), 
                       rad_range=c(1, 1), noise=F, r_noise=3, n_noise=15){
  
  cluster_data = lapply(1:nrow(centers), function(i){
    center = centers[i, ] # the center of the current cluster
    radius = runif(1, rad_range[1], rad_range[2]) # the random radius
    cluster = rpoisball_unit(n_cls[i], d)*radius+
      matrix(rep(center, n_cls[i]), nrow = n_cls[i], byrow=T) # generate a cluster
    return(cluster)
  })
  
  cluster_data = do.call(rbind, cluster_data)
  
  if(noise){
    center_noise = colMeans(centers)
    noise_data = rpoisball_unit(n_noise, d)*r_noise+center_noise
    cluster_data = rbind(cluster_data, noise_data)
  }
  
  return(cluster_data)
}

# simulate a data set with Gaussian clusters and noise (no outliers)
# n_cls: the number of each cluster
# d: dimension
# center: the center of each cluster by row
# rad_range: the upper & lower bound of random radius
# noise: whether there is noise
# r_noise: the radius of noise cluster
# n_noise: the number of noise
Gaussian_cls = function(n_cls=c(100,100), d=2, 
                       centers=rbind(c(0, 0), c(0, 3)), sigma=0.35,
                       rad_range=c(1, 1), noise=F, r_noise=3, n_noise=15){
  Sigma = diag(sigma, d)
  cluster_data = lapply(1:nrow(centers), function(i){
    center = centers[i, ] # the center of the current cluster
    radius = runif(1, rad_range[1], rad_range[2]) # the random radius
    cluster = mvrnorm(n_cls[i], center, Sigma*(radius^2)) # generate a cluster
    return(cluster)
  })
  
  cluster_data = do.call(rbind, cluster_data)
  
  if(noise){
    center_noise = colMeans(centers)
    noise_data = rpoisball_unit(n_noise, d)*r_noise+center_noise
    cluster_data = rbind(cluster_data, noise_data)
  }
  
  return(cluster_data)
}

# realisation the mixture of the Matern Cluster process and the Thomas process inside a (hyper) unit square
# d = dimension
# kappa1: Intensity of the Poisson process of cluster centres of the Matérn Cluster process
# r: Radius of the clusters. A single positive number.
# mu1: Mean number of points per cluster (a single positive number)
# expand1: Window expansion distance of the Matérn Cluster process for parents

# kappa2: Intensity of the Poisson process of cluster centres of the Thomas process
# scale: Cluster size. Standard deviation of random displacement (along each coordinate axis) of a point from its cluster centre. A single positive number.
# mu2: Mean number of points per cluster (a single positive number)
# slen: the length of the square window.
# expand2: Window expansion distance of the Thomas process for parents
# nsim: Number of simulated realisations to be generated.
# kappa_O: the intensity of outliers

Uni.Gau_cls = function(d, kappa1, r, mu1, expand1=0, kappa2, scale, mu2, expand2=0, slen=1, kappa_O){
  n1 = 0
  n2 = 0
  n0 = 0
  while((n1+n2)<80 | n0 == 0){
    # Generate the Matérn Cluster process
    np1 = rpois(1,kappa1*(slen+2*expand1)^d) # simulate the number of parents
    parents1 = matrix(runif(np1*d,min=-expand1,max=slen+expand1),ncol=d,byrow=T) # simulate the location of parents
    
    cls1 = apply(parents1,1,function(parent){ # simulate each children cluster around parents
      n_child = rpois(1, mu1)  # simulate the number of cluster
      if(n_child!=0){
        child = rpoisball.unit(n_child,d)*r + matrix(rep(parent,n_child),ncol=d,byrow=T)
      } else {child = NULL}
      return(child)
    },simplify=F)
    
    if(length(cls1)>0){
      cls1 = do.call(rbind, cls1)
      index = sapply(1:nrow(cls1),function(x){
        if(all(cls1[x,]>0) & all(cls1[x,]<slen)) return(x)})
      index = unlist(index)
      cls1 = cls1[index,]
    }
    n1 = length(cls1)/d
    
    # Generate the Thomas process
    np2 = rpois(1,kappa2*(slen+2*expand2)^d) # the number of parents
    parents2 = matrix(runif(np2*d,min=-expand2,max=slen+expand2),ncol=d,byrow=T) # the location of parents
    
    cls2 = apply(parents2,1,function(parent){  # simulate each children cluster
      n_child = rpois(1, mu2)
      if(n_child!=0){
        child = mvrnorm(n_child,parent,diag(d)*scale)
      } else {child = NULL}
      return(child)
    },simplify=F)
    
    if(length(cls2)>0){
      cls2 = do.call(rbind, cls2)
      index = sapply(1:nrow(cls2),function(x){
        if(all(cls2[x,]>0) & all(cls2[x,]<slen)) return(x)
      })
      index = unlist(index)
      cls2 = cls2[index,]
    }
    n2 = length(cls2)/d
    
    index = NULL
    np3 = 0
    # Generate outliers
    while(length(index)==0 | np3 <= 1){
      np3 = rpois(1,kappa_O*slen^d) # simulate the number of outliers
      outlier = matrix(runif(np3*d,min=0,max=slen),ncol=d,byrow=T)
      
      index = sapply(1:nrow(outlier),function(x){
        dist.mat1 = t(t(parents1)-outlier[x,])
        dist1 = sqrt(apply(dist.mat1^2,1,sum))
        
        dist.mat2 = t(t(parents2)-outlier[x,])
        dist2 = sqrt(apply(dist.mat2^2,1,sum))
        
        if(all(dist1>=(2*r)) & all(dist2>=(3.33*sqrt(scale)))) return(x)
      })
      index = unlist(index)
    }
    
    cls3 = outlier[-index,] # extract "outliers" that are not too far from parents
    outlier = outlier[index,]
    n3 = length(cls3)/d 
    n0 = length(outlier)/d
  }
    return(c(Matérn_parents = list(parents1), Matérn_children = list(cls1), 
             Thomas_parents = list(parents2), Thomas_children = list(cls2),
             noise = list(cls3),
             Outlier = list(outlier), 
             num = list(c(n1+n2+n3,n0))))
}