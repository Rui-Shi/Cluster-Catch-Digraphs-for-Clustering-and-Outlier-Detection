# x is the data set feeded
# C is a penalty parameter controls the value of the radius
svdd_sphere = function(x,C=0.05){
  library(reticulate)
  py_run_string("import os")
  py_run_string("os.chdir('D:/code_working_folder/KS-MCGs/simulations/Uniform_MeanSlippage/SVDD-Python-master')")
  py_run_string("import sys")
  py_run_string("import numpy as np")
  py_run_string("from src.BaseSVDD import BaseSVDD")
  # svdd object using polynomial kernel with degree 1, which has rigorous sphere boundary
  py_run_string("svdd = BaseSVDD(C=0.1, kernel='poly',degree=1, display='off')")
  py$data = x
  py_run_string("svdd.fit(data)")
  py_run_string("radius = svdd.radius")
  py_run_string("center = svdd.center")
  return(list(center = py$center, radius = py$radius))
}

svdd_sphere

