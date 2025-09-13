source("G:/code_working_folder/general functions/Uni-Gau_cls.R")

# plot some simulations of Matern Cluster process
ratio = c(1.000000, 1.024079, 1.068926, 1.163903, 1.291596)

d = 2
iteN = 1000
quant = 0.99

# simulation settings
kappa1 = 6
mu1 = 36.4*ratio[1]
expand1 = 0
r = 0.1
kappa2 = 0
scale = 0.005
mu2 = 36.4*ratio[1]
expand2 = 0
slen = 1
kappa_O = 20
# simulate clusters of random sizes and positions
set.seed(1234)
data.list = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  return(data_simu)
})

# 1
par(cex=1)
result = data.list[[2]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 2
result = data.list[[20]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 3
result = data.list[[200]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 4
result = data.list[[1000]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)



# plot some simulations of Thomas Cluster process
ratio = c(1.000000, 1.024079, 1.068926, 1.163903, 1.291596)

d = 2
iteN = 1000
quant = 0.99

# simulation settings
kappa1 = 6
mu1 = 36.4*ratio[1]
expand1 = 0
r = 0.1
kappa2 = 0
scale = 0.005
mu2 = 36.4*ratio[1]
expand2 = 0
slen = 1
kappa_O = 20
# simulate clusters of random sizes and positions
set.seed(1234)
data.list = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  return(data_simu)
})

# 1
par(cex=1)
result = data.list[[2]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 2
result = data.list[[20]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 3
result = data.list[[200]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 4
result = data.list[[1000]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)




# plot some simulations of Thomas Cluster process
ratio = c(1.000000, 1.064383, 1.193162, 1.604715, 2.859147)

d = 2
iteN = 1000
quant = 0.99

# simulation settings
kappa1 = 0
mu1 = 37.1*ratio[1]
expand1 = 0
r = 0.1
kappa2 = 6
scale = 0.005
mu2 = 37.1*ratio[1]
expand2 = 0
slen = 1
kappa_O = 20
# simulate clusters of random sizes and positions
set.seed(1234)
data.list = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  return(data_simu)
})

# 1
par(cex=1)
result = data.list[[3]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 2
result = data.list[[30]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 3
result = data.list[[300]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 4
result = data.list[[1000]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)



# plot some simulations of Thomas Cluster process
ratio = c(1.000000, 1.024079, 1.068926, 1.163903, 1.291596)

d = 2
iteN = 1000
quant = 0.99

# simulation settings
kappa1 = 6
mu1 = 36.4*ratio[1]
expand1 = 0
r = 0.1
kappa2 = 0
scale = 0.005
mu2 = 36.4*ratio[1]
expand2 = 0
slen = 1
kappa_O = 20
# simulate clusters of random sizes and positions
set.seed(1234)
data.list = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  return(data_simu)
})

# 1
par(cex=1)
result = data.list[[2]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 2
result = data.list[[20]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 3
result = data.list[[200]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 4
result = data.list[[1000]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)


# plot some simulations of Thomas Cluster process
ratio = c(1.000000, 1.064383, 1.193162, 1.604715, 2.859147)

d = 2
iteN = 1000
quant = 0.99

# simulation settings
kappa1 = 0
mu1 = 37.1*ratio[1]
expand1 = 0
r = 0.1
kappa2 = 6
scale = 0.005
mu2 = 37.1*ratio[1]
expand2 = 0
slen = 1
kappa_O = 20
# simulate clusters of random sizes and positions
set.seed(1234)
data.list = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  return(data_simu)
})

# 1
par(cex=1)
result = data.list[[3]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 2
result = data.list[[30]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 3
result = data.list[[300]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 4
result = data.list[[1000]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)



# plot some simulations of Thomas Cluster process
ratio = c(1.000000, 1.051859, 1.137064, 1.338852, 1.772706)

d = 2
iteN = 1000
quant = 0.99

# simulation settings
kappa1 = 3
mu1 = 36.7*ratio[1]
expand1 = 0
r = 0.1
kappa2 = 3
scale = 0.005
mu2 = 36.7*ratio[1]
expand2 = 0
slen = 1
kappa_O = 20
# simulate clusters of random sizes and positions
set.seed(1234)
data.list = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  return(data_simu)
})

# 1
par(cex=1)
result = data.list[[2]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 2
result = data.list[[40]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 3
result = data.list[[200]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)

# 4
result = data.list[[1000]]
plot(rbind(result$Matérn_children,result$Thomas_children),asp=1,
     ylim=c(0,1),xlim=c(0,1),pch=20,ylab="",xlab="",axes=FALSE,frame.plot = T)
points(result$Outlier,col="red",pch=20)
points(result$Matérn_parents,col="green",pch=19)
points(result$Thomas_parents,col="blue",pch=19)