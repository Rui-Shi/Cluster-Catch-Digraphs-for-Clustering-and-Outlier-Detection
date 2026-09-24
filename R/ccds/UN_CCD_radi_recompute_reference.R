# Verbatim copy of nnccd.radi() from R/ccds/UN_CCD.R as it stood before the
# incremental nearest-neighbour update (2026-09-23), renamed nnccd.radi.recompute.
# Kept only as the reference for revision_experiments/tr1/92_validate_incremental_radi.R.
# Nothing in the detectors sources this file.

nnccd.radi.recompute <- function(dx, quantile="lower", method="ascend", low.num, quant, simul=NULL, niter, scores=F){
  
  ddx <- as.matrix(dist(dx)) # the distance matrix
  n <- nrow(dx)
  d <- ncol(dx)
  R <- rep(0,n)
  
  if(quantile=="lower"){
    if(!is.null(simul)) {NN.envelop <- list(average=simul$average[1:n],median=simul$median[1:n])} 
    else {NN.envelop <- NNDest.simpois.lower.quant(n, d, quant, niter)}
    if(!scores){
      for(i in 1:n){
        if(method == "ascend"){
          o.d <- order(ddx[i,]) # the descending distance order for i_th object
          for(j in low.num:n){
            r <- ddx[i,o.d[j]]
            NN.dist.obs <- NNDest.dist.f(ddx[o.d[2:j],o.d[2:j]],r) # the average NN distance of within a covering ball, the center point is dropped
            
            # check the values, if accepted, set the R[i] as the radius
            lower.bound.ave = NN.envelop$average[j-1]
            lower.bound.med = NN.envelop$median[j-1]
            # if(NN.dist.obs$averge<lower.bound.ave | NN.dist.obs$median<lower.bound.med){
            #   R[i] = ddx[i,o.d[j-1]]
            #   break
            # }
            if(NN.dist.obs$averge<lower.bound.ave | NN.dist.obs$median<lower.bound.med){
              if(j == low.num) R[i] = 0
              else  R[i] = ddx[i,o.d[j-1]]
              break
            }
          }
        }
        if(method=="descend"){
          o.d <- order(ddx[i,], decreasing=T) # the descending distance order for i_th object
          for(j in 1:(n-low.num)){
            r <- ddx[i,o.d[j]]
            NN.dist.obs <- NNDest.dist.f(ddx[o.d[j:(n-1)],o.d[j:(n-1)]],r) # the average NN distance of within a covering ball, the center point is dropped
            
            # check the values, if accepted, set the R[i] as the radius
            lower.bound.ave = rev(NN.envelop$average)[j+2]
            lower.bound.med = rev(NN.envelop$median)[j+2]
            if(NN.dist.obs$averge>lower.bound.ave & NN.dist.obs$median>lower.bound.med){
              R[i] = r
              break
            }
          }
        }
      }
    } else {
      for(i in 1:n){
        if(method == "ascend"){
          o.d <- order(ddx[i,]) # the descending distance order for i_th object
          for(j in low.num:n){
            r <- ddx[i,o.d[j]]
            NN.dist.obs <- NNDest.dist.f(ddx[o.d[2:j],o.d[2:j]],r) # the average NN distance of within a covering ball, the center point is dropped
            
            # check the values, if accepted, set the R[i] as the radius
            lower.bound.ave = NN.envelop$average[j-1]
            lower.bound.med = NN.envelop$median[j-1]
            # if(NN.dist.obs$averge<lower.bound.ave | NN.dist.obs$median<lower.bound.med){
            #   R[i] = ddx[i,o.d[j-1]]
            #   break
            # }
            if(NN.dist.obs$averge<lower.bound.ave | NN.dist.obs$median<lower.bound.med){
              if(j == low.num) R[i] = 0
              else  R[i] = ddx[i,o.d[j-1]]
              break
            }
          }
        }
        if(method=="descend"){
          o.d <- order(ddx[i,], decreasing=T) # the descending distance order for i_th object
          for(j in 1:(n-low.num)){
            r <- ddx[i,o.d[j]]
            NN.dist.obs <- NNDest.dist.f(ddx[o.d[j:(n-1)],o.d[j:(n-1)]],r) # the average NN distance of within a covering ball, the center point is dropped
            
            # check the values, if accepted, set the R[i] as the radius
            lower.bound.ave = rev(NN.envelop$average)[j+2]
            lower.bound.med = rev(NN.envelop$median)[j+2]
            if(NN.dist.obs$averge>lower.bound.ave & NN.dist.obs$median>lower.bound.med){
              R[i] = r
              break
            }
          }
        }
        if(R[i]==0){R[i]=sort(ddx[i,])[2]} # avoid 0 radius (necessary for outlyingness scores!)
      }
    }
  }
  return(list(R=R,KS=NULL))
}
