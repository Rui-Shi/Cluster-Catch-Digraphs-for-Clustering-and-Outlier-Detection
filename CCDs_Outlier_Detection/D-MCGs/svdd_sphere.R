# x is the data set feeded
# C is a penalty parameter controls the value of the radius
svdd_sphere = function(x,C=0.05){
  library(reticulate)
  py_run_string("import os")
  py_run_string("os.chdir('/mmfs1/home/rzs0112/code_working_folder/KS-MCGs/SVDD-Python-master')")
  py_run_string("import sys")
  py_run_string("import numpy as np")
  py_run_string("from src.BaseSVDD import BaseSVDD")
  # svdd object using polynomial kernel with degree 1, which has rigorous sphere boundary
  py$P = C
  py_run_string("svdd = BaseSVDD(C=P, kernel='poly',degree=1, display='off')")
  py$data = x
  py_run_string("svdd.fit(data)")
  py_run_string("radius = svdd.radius")
  py_run_string("center = svdd.center")
  return(list(center = py$center, radius = py$radius))
}

