#####
# m: record the maximal density parameter m such that the KS-CCD is connected.
# step: the value m decreases each time
# m.intial: the starting value of m
source(here::here("R/ccds/ccd_ks_NEW.R"))
source(here::here("R/ccds/functions.R"))
connected.ksccd.m = function(t){
  source(here::here("R/ccds/ccd_ks_NEW.R"))
  source(here::here("R/ccds/functions.R"))
  #source(here::here("R/ccds/ccd_ks_NEW.R"))
  #source(here::here("R/ccds/functions.R"))
  #source(here::here("R/ccds/ccd_ks_NEW.R"))
  #source(here::here("R/ccds/functions.R"))
  step = 100
  m.intial = 100
  member = c(1,2)
  m = m.intial
  while(length(unique(member))>1){
    m = m - step
    member = ksccd.connected(t,m,sequential=FALSE,alpha=0.05)$member
  }
  
  step = step/10 #smaller steps
  while(length(unique(member))==1){
    m = m + step
    member = ksccd.connected(t,m,sequential=FALSE,alpha=0.05)$member
  }
  
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
  return(m)
}