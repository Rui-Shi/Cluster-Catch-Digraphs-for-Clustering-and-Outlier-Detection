library(igraph)
#####
# m: record the maximal density parameter m such that the KS-CCD is connected.
# step: the value m decreases each time
# m.intial: the starting value of m
connected.ksccd.m.old = function(t){
  source(here::here("R/ccds/ccd_ks_NEW.R"))
  source(here::here("R/ccds/ccdfunctions.R"))
  #source(here::here("R/ccds/ccd_ks_NEW.R"))
  #source(here::here("R/ccds/functions.R"))
  #source(here::here("R/ccds/ccd_ks_NEW.R"))
  #source(here::here("R/ccds/functions.R"))
  if(is.null(nrow(t))){
    
    m = 1
    
  } else {
    step = 100
    m.intial = 100
    member = c(1,2)
    m = m.intial
    i = -1
    while(m/step<=10){
      step = step/10 #smaller steps
      while(length(unique(member))>1){
        m = round(m - step,i)
        member = ksccd.connected(t,m,sequential=FALSE,alpha=0.05)$member
      }
      
      step = step/10 #smaller steps
      while(length(unique(member))==1){
        m = round(m + step,i+1)
        member = ksccd.connected(t,m,sequential=FALSE,alpha=0.05)$member
      }
      i = i + 2
      m = m - step #reverse the last step
    }
  }
  return(m)
}





#####
# m: record the maximal density parameter m such that the KS-CCD is connected.
# step: the value m decreases each time
# m.intial: the starting value of m
connected.ksccd.m = function(t){
   source(here::here("R/ccds/ccd_ks_NEW.R"))
  source(here::here("R/ccds/ccdfunctions.R"))
  # source(here::here("R/ccds/ccd_ks_NEW.R"))
  # source(here::here("R/ccds/ccdfunctions.R"))
  #   source(here::here("R/ccds/ccd_ks_NEW.R"))
  #   source(here::here("R/ccds/ccdfunctions.R"))
  if(is.null(nrow(t))){
    m = 1
  } else {
    m.intial = 1
    m = m.intial
    member = ksccd.connected(t,m,sequential=FALSE,alpha=0.05)$member
    while(length(unique(member))==1){
      m = m*10
      member = ksccd.connected(t,m,sequential=FALSE,alpha=0.05)$member
    }
    
    step = m
    
    while(m/step<=10 | step>0.01){
      
      step = step/10 #smaller steps
      while(length(unique(member))>1){
        m = m - step
        member = ksccd.connected(t,m,sequential=FALSE,alpha=0.05)$member
      }
      
      step = step/10 #smaller steps
      while(length(unique(member))==1){
        m = m + step
        member = ksccd.connected(t,m,sequential=FALSE,alpha=0.05)$member
      }
      m = m - step #reverse the last step
      if(m/step>=1e+12) break
    }
  }
  return(m)
}

