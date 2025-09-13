library(e1071)
library(quadprog)
library(nloptr)


objfun = function(x, data){ # objective function L
  data_alpha = diag(x)%*%data
  L = sum(x*diag(data%*%t(data))) - sum(data_alpha%*%t(data_alpha))
  return(-L)
}

center_list = NULL
for(i in 1:100){
  datax = data.list[[i]][1:48,]
  n = dim(datax)[1] # length of the data set
  C = 0.1 # penality parameter
  x0 = rep(C/2,n) # vector with staring values for alphas
  
  res = constrOptim(theta = x0, ui = rep(1,n),  f= objfun, data = datax, method = "L-BFGS-B", lower=rep(0,n), upper=rep(C,n)) 
  print(res$par)
  center = colSums(res$par*datax)/sum(res$par)
  center_list = rbind(center_list,center)
}

c(1,2)*matrix(c(1,2,3,4),nrow=2)




objfun = function(x){ # objective function L
  data_alpha = diag(x)%*%data
  L = sum(x*diag(data%*%t(data))) - sum(data_alpha%*%t(data_alpha))
  return(-L)
}
C = 0.1 # penality parameter
data = data.list[[6]]
n = dim(data)[1] # length of the data set
library(nloptr)
# Define the constraint function
confun <- function(x) {
  return(sum(x)-1)
}

# Set the optimization options
local_opts <- list( "algorithm" = "NLOPT_GN_ISRES", "xtol_rel" = 1.0e-4 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-4,
              "maxeval"= 20000000,
              "local_opts" = local_opts,
              "print_level" = 0 )

# Run the optimization
result <- nloptr(x0 = c(rep(C,n)),
                 eval_f = objfun,
                 lb = rep(0,n),
                 ub = rep(C,n),
                 eval_g_eq = confun,
                 opts = opts,
                 )

result$solution
sum(result$solution)

colSums(result$solution*data)
#This code uses the nloptr package in R to solve the problem. The objfun variable represents the objective function. The confun variable represents the constraint function. The lb and ub variables represent the lower and upper bounds for the variables. The g variable represents the constraints. The opts variable represents the optimization options.

#You can modify this code to fit your specific problem by changing the values of objfun, confun, lb, ub, g, and opts.

#I hope this helps!


M = matrix(c(1,2,3,4),byrow=T,nrow=2)
M%*%t(M)
c(1,2)*diag(M)