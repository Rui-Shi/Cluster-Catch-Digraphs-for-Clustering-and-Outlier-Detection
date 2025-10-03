source("/media/rui/exNVME/code_working_folder/general functions/Uni-Gau_cls.R")
d = 2
kappa1 = 0
mu1 = 33.7
expand1 = 0
r = 0.1
kappa2 = 6
scale = 0.005
mu2 = 33.7
expand2 = 0
slen = 1
kappa_O = 20
iteN = 1000

# simulate clusters of random sizes and positions
set.seed(1234)
data.listNum2d = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  cls1 = data_simu$Matérn_children
  cls2 = data_simu$Thomas_children
  outlier = data_simu$Outlier
  return(c(data = list(rbind(cls1, cls2, outlier)), num = list(data_simu$num)))
})

data.list2d = lapply(1:iteN, function(x){
  return(data.listNum2d[[x]]$data)
})

data.num2d = lapply(1:iteN, function(x){
  return(data.listNum2d[[x]]$num)
})
data.num2d = do.call(rbind, data.num2d)
mean2d = apply(data.num2d, 2, mean)
print(mean2d)


d = 3
set.seed(123)
data.listNum3d = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  cls1 = data_simu$Matérn_children
  cls2 = data_simu$Thomas_children
  outlier = data_simu$Outlier
  return(c(data = list(rbind(cls1, cls2, outlier)), num = list(data_simu$num)))
})

data.list3d = lapply(1:iteN, function(x){
  return(data.listNum3d[[x]]$data)
})

data.num3d = lapply(1:iteN, function(x){
  return(data.listNum3d[[x]]$num)
})
data.num3d = do.call(rbind, data.num3d)
mean3d = apply(data.num3d, 2, mean)
print(mean3d)


d = 5
set.seed(12)
data.listNum5d = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  cls1 = data_simu$Matérn_children
  cls2 = data_simu$Thomas_children
  outlier = data_simu$Outlier
  return(c(data = list(rbind(cls1, cls2, outlier)), num = list(data_simu$num)))
})

data.list5d = lapply(1:iteN, function(x){
  return(data.listNum5d[[x]]$data)
})

data.num5d = lapply(1:iteN, function(x){
  return(data.listNum5d[[x]]$num)
})
data.num5d = do.call(rbind, data.num5d)
mean5d = apply(data.num5d, 2, mean)
print(mean5d)


d = 10
set.seed(12)
data.listNum10d = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  cls1 = data_simu$Matérn_children
  cls2 = data_simu$Thomas_children
  outlier = data_simu$Outlier
  return(c(data = list(rbind(cls1, cls2, outlier)), num = list(data_simu$num)))
})

data.list10d = lapply(1:iteN, function(x){
  return(data.listNum10d[[x]]$data)
})

data.num10d = lapply(1:iteN, function(x){
  return(data.listNum10d[[x]]$num)
})
data.num10d = do.call(rbind, data.num10d)
mean10d = apply(data.num10d, 2, mean)
print(mean10d)


d = 20
set.seed(12)
data.listNum20d = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  cls1 = data_simu$Matérn_children
  cls2 = data_simu$Thomas_children
  outlier = data_simu$Outlier
  return(c(data = list(rbind(cls1, cls2, outlier)), num = list(data_simu$num)))
})

data.list20d = lapply(1:iteN, function(x){
  return(data.listNum20d[[x]]$data)
})

data.num20d = lapply(1:iteN, function(x){
  return(data.listNum20d[[x]]$num)
})
data.num20d = do.call(rbind, data.num20d)
mean20d = apply(data.num20d, 2, mean)
print(mean20d)

mu1 = mu1*20
mu2 = mu2*20
d = 50
set.seed(14)
data.listNum50d = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  cls1 = data_simu$Matérn_children
  cls2 = data_simu$Thomas_children
  outlier = data_simu$Outlier
  return(c(data = list(rbind(cls1, cls2, outlier)), num = list(data_simu$num)))
})

data.list50d = lapply(1:iteN, function(x){
  return(data.listNum50d[[x]]$data)
})

data.num50d = lapply(1:iteN, function(x){
  return(data.listNum50d[[x]]$num)
})
data.num50d = do.call(rbind, data.num50d)
mean50d = apply(data.num50d, 2, mean)
print(mean50d)

mu1 = mu1*20
mu2 = mu2*20
d = 100
set.seed(14)
data.listNum100d = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  cls1 = data_simu$Matérn_children
  cls2 = data_simu$Thomas_children
  outlier = data_simu$Outlier
  return(c(data = list(rbind(cls1, cls2, outlier)), num = list(data_simu$num)))
})

data.list100d = lapply(1:iteN, function(x){
  return(data.listNum100d[[x]]$data)
})

data.num100d = lapply(1:iteN, function(x){
  return(data.listNum100d[[x]]$num)
})
data.num100d = do.call(rbind, data.num100d)
mean100d = apply(data.num100d, 2, mean)
print(mean100d)


ratio = c(33.7, 36.12637, 40.89483*200/193, 50.19943*200/182, 65.64452*200/145, 467.29827*200/182, 6419.09588*200/175)
print(ratio)




# after adjust
d = 2
kappa1 = 0
mu1 = 33.7
expand1 = 0
r = 0.1
kappa2 = 6
scale = 0.005
mu2 = 33.7
expand2 = 0
slen = 1
kappa_O = 20
iteN = 1000

# simulate clusters of random sizes and positions
set.seed(1234)
data.listNum2d = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  cls1 = data_simu$Matérn_children
  cls2 = data_simu$Thomas_children
  outlier = data_simu$Outlier
  return(c(data = list(rbind(cls1, cls2, outlier)), num = list(data_simu$num)))
})

data.list2d = lapply(1:iteN, function(x){
  return(data.listNum2d[[x]]$data)
})

data.num2d = lapply(1:iteN, function(x){
  return(data.listNum2d[[x]]$num)
})
data.num2d = do.call(rbind, data.num2d)
mean2d = apply(data.num2d, 2, mean)
print(mean2d)


mu1 = ratio[2]
mu2 = ratio[2]
d = 3
set.seed(123)
data.listNum3d = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  cls1 = data_simu$Matérn_children
  cls2 = data_simu$Thomas_children
  outlier = data_simu$Outlier
  return(c(data = list(rbind(cls1, cls2, outlier)), num = list(data_simu$num)))
})

data.list3d = lapply(1:iteN, function(x){
  return(data.listNum3d[[x]]$data)
})

data.num3d = lapply(1:iteN, function(x){
  return(data.listNum3d[[x]]$num)
})
data.num3d = do.call(rbind, data.num3d)
mean3d = apply(data.num3d, 2, mean)
print(mean3d)


mu1 = ratio[3]
mu2 = ratio[3]
d = 5
set.seed(12)
data.listNum5d = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  cls1 = data_simu$Matérn_children
  cls2 = data_simu$Thomas_children
  outlier = data_simu$Outlier
  return(c(data = list(rbind(cls1, cls2, outlier)), num = list(data_simu$num)))
})

data.list5d = lapply(1:iteN, function(x){
  return(data.listNum5d[[x]]$data)
})

data.num5d = lapply(1:iteN, function(x){
  return(data.listNum5d[[x]]$num)
})
data.num5d = do.call(rbind, data.num5d)
mean5d = apply(data.num5d, 2, mean)
print(mean5d)



mu1 = ratio[4]
mu2 = ratio[4]
d = 10
set.seed(12)
data.listNum10d = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  cls1 = data_simu$Matérn_children
  cls2 = data_simu$Thomas_children
  outlier = data_simu$Outlier
  return(c(data = list(rbind(cls1, cls2, outlier)), num = list(data_simu$num)))
})

data.list10d = lapply(1:iteN, function(x){
  return(data.listNum10d[[x]]$data)
})

data.num10d = lapply(1:iteN, function(x){
  return(data.listNum10d[[x]]$num)
})
data.num10d = do.call(rbind, data.num10d)
mean10d = apply(data.num10d, 2, mean)
print(mean10d)


mu1 = ratio[5]
mu2 = ratio[5]
d = 20
set.seed(12)
data.listNum20d = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  cls1 = data_simu$Matérn_children
  cls2 = data_simu$Thomas_children
  outlier = data_simu$Outlier
  return(c(data = list(rbind(cls1, cls2, outlier)), num = list(data_simu$num)))
})

data.list20d = lapply(1:iteN, function(x){
  return(data.listNum20d[[x]]$data)
})

data.num20d = lapply(1:iteN, function(x){
  return(data.listNum20d[[x]]$num)
})
data.num20d = do.call(rbind, data.num20d)
mean20d = apply(data.num20d, 2, mean)
print(mean20d)


mu1 = ratio[6]
mu2 = ratio[6]
d = 50
set.seed(12)
data.listNum50d = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  cls1 = data_simu$Matérn_children
  cls2 = data_simu$Thomas_children
  outlier = data_simu$Outlier
  return(c(data = list(rbind(cls1, cls2, outlier)), num = list(data_simu$num)))
})

data.list50d = lapply(1:iteN, function(x){
  return(data.listNum50d[[x]]$data)
})

data.num50d = lapply(1:iteN, function(x){
  return(data.listNum50d[[x]]$num)
})
data.num50d = do.call(rbind, data.num50d)
mean50d = apply(data.num50d, 2, mean)
print(mean50d)


mu1 = ratio[7]
mu2 = ratio[7]
d = 100
set.seed(14)
data.listNum100d = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  cls1 = data_simu$Matérn_children
  cls2 = data_simu$Thomas_children
  outlier = data_simu$Outlier
  return(c(data = list(rbind(cls1, cls2, outlier)), num = list(data_simu$num)))
})

data.list100d = lapply(1:iteN, function(x){
  return(data.listNum100d[[x]]$data)
})

data.num100d = lapply(1:iteN, function(x){
  return(data.listNum100d[[x]]$num)
})
data.num100d = do.call(rbind, data.num100d)
mean100d = apply(data.num100d, 2, mean)
print(mean100d)
