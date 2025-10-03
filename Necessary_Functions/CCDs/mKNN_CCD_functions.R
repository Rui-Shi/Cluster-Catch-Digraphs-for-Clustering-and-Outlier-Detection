library(igraph)
#####
# m: record the maximal density parameter m such that the KS-CCD is connected.
# step: the value m decreases each time
# m.intial: the starting value of m
connected.ksccd.m.old = function(t){
  source("/mmfs1/home/rzs0112/code_working_folder/ccds/ccd_ks_NEW.R")
  source("/mmfs1/home/rzs0112/code_working_folder/ccds/ccdfunctions.R")
  #source("C:/Users/shiru/OneDrive - Auburn University/Research Outliers Detection/code/ccds/ccd_ks_NEW.R")
  #source("C:/Users/shiru/OneDrive - Auburn University/Research Outliers Detection/code/ccds/functions.R")
  #source("G:/OneDrive - Auburn University/Research Outliers Detection/code/ccds/ccd_ks_NEW.R")
  #source("G:/OneDrive - Auburn University/Research Outliers Detection/code/ccds/functions.R")
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
   source("/mmfs1/home/rzs0112/code_working_folder/ccds/ccd_ks_NEW.R")
  source("/mmfs1/home/rzs0112/code_working_folder/ccds/ccdfunctions.R")
  # source("C:/Users/shiru/OneDrive - Auburn University/Research Outliers Detection/code/ccds/ccd_ks_NEW.R")
  # source("C:/Users/shiru/OneDrive - Auburn University/Research Outliers Detection/code/ccds/ccdfunctions.R")
  #   source("/media/rui/exNVME/code_working_folder/ccds/ccd_ks_NEW.R")
  #   source("/media/rui/exNVME/code_working_folder/ccds/ccdfunctions.R")
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

