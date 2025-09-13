library(plotrix)
library(dplyr)
library(MASS)
#Data Simulation
# We have 3 clusters centered at (3,13),(16,15),(10,10)
# We also has 3 outliers: (10,17),(5,5), (17,5)
set.seed(1)
mean1 = c(3,13)
mean2 = c(16,15)
mean3 = c(10,10)
sigma1 = matrix(c(1,0,0,1),nrow=2)
sigma2 = matrix(c(2,1,1,2),nrow=2)
sigma3 = matrix(c(3,2,2,3),nrow=2)
c1 = mvrnorm(100,mean1,sigma1)
c2 = mvrnorm(120,mean2,sigma2)
c3 = mvrnorm(150,mean3,sigma3)
o1 = c(10,17)
o2 = c(5,5)
o3 = c(17,5)

par(cex=0.7)
data = rbind(c1,c2,c3,o1,o2,o3) # combining simulated data
plot(data[1:100,1],data[1:100,2],xlab = "X", ylab = "Y",pch=20,ylim = c(2,20),xlim=c(-3,22),col="red")
points(data[101:220,1],data[101:220,2],xlab = "X", ylab = "Y",pch=20,col="darkgreen")
points(data[221:370,1],data[221:370,2],xlab = "X", ylab = "Y",pch=20,col="blue")
points(data[371:373,1],data[371:373,2],xlab = "X", ylab = "Y",pch=3)

legend("topleft",                                  
         legend=c("Cluster1","Cluster2","Cluster3","Outliers"),        
         col=c("red","darkgreen","blue","black"),                 
         pch=c(20,20,20,3),cex=1.5)                                         


plot(data[1:100,1],data[1:100,2],xlab = "X", ylab = "Y",pch=20,ylim = c(2,20),xlim=c(-3,22),col="red")
points(data[101:220,1],data[101:220,2],xlab = "X", ylab = "Y",pch=20,col="darkgreen")
points(data[221:370,1],data[221:370,2],xlab = "X", ylab = "Y",pch=20,col="blue")
points(data[371:373,1],data[371:373,2],xlab = "X", ylab = "Y",pch=3)

legend("topleft",                                  
       legend=c("Cluster1","Cluster2","Cluster3","Outliers"),        
       col=c("red","darkgreen","blue","black"),                 
       pch=c(20,20,20,3),cex=1.5)


############
##K-S CCDs##
############

delta = 10 #Set density parameter to be 10
d = 2 #dimension
dist.mat = NULL #distance matrix of the whole data set
len = length(data[,1]) #length of the data set

cover.num = NULL # record the number of points covered by each covering ball
cover = NULL # a list stores the objects covered by each object
r.opt = NULL # optimal radius for each objects

for(a in 1:len){
  dist = NULL
  for(b in 1:len){
    dist.temp = sqrt((data[a,1]-data[b,1])^2+(data[a,2]-data[b,2])^2) #dist.temp: the distance between point a and b
    dist = c(dist,dist.temp) # dist is a vector representing the distances of a with other points in the data set 
  }
  dist.mat = rbind(dist.mat,dist)  # determine the distance between the current point and all others.
}
rownames(dist.mat)=NULL
colnames(dist.mat)=NULL

#determine the radius of the covering ball for each points by K-S statistic
ks.seq = NULL # record K-S statistic for all the objects
for(a in 1:len){
  r.seq = sort(dist.mat[a,]) 
  ks.temp = NULL # record k-s statistics
  for(i in 1:len){
    r.temp = r.seq[i]
    ks.ttemp = length(which(dist.mat[a,]<r.temp))-delta*(r.temp^d)
    ks.temp=c(ks.temp,ks.ttemp)
  }
  ks.seq = c(ks.seq, max(ks.temp))
  index.max = which.max(ks.temp)
  cover.num = c(cover.num,length(which(dist.mat[a,]<r.seq[index.max])))
  cover = c(cover,list(which(dist.mat[a,]<r.seq[index.max])))
  r.opt = c(r.opt,r.seq[index.max])
  draw.circle(data[a,1],data[a,2],r.seq[index.max],nv=1000,lwd=0.01)
}
text(17,3,expression(paste(delta, " = 5")),cex=1.5)
names(r.opt)=NULL


# Finding an approximate minimal dominating set (algoritm 2 in Parameter Free Clustering with CCDs)
plot(data[1:100,1],data[1:100,2],xlab = "X", ylab = "Y",pch=20,ylim = c(2,20),xlim=c(-3,22),col="red",main="The Covering Balls of An Approx MDS (alg 2)")
points(data[101:220,1],data[101:220,2],xlab = "X", ylab = "Y",pch=20,col="darkgreen")
points(data[221:370,1],data[221:370,2],xlab = "X", ylab = "Y",pch=20,col="blue")
points(data[371:373,1],data[371:373,2],xlab = "X", ylab = "Y",pch=3)
text(17,3,expression(paste(delta, " = 5")),cex=1.5)
legend("topleft",                                  
       legend=c("Cluster1","Cluster2","Cluster3","Outliers"),        
       col=c("red","darkgreen","blue","black"),                 
       pch=c(20,20,20,3),cex=1.5)  

MDS = NULL # store the elements in approximate MDS
r.opt = r.opt # optimal radius for MDS
dist.matM = dist.mat #the distance matrix used for finding the MDS
len0 = length(dist.matM[1,])
len1 = length(dist.matM[,1])
cover.MDS = NULL


while(len0>0){
  cover.num1 = NULL #re-calculate the covering number of each balls in each step
  for(i in 1:len1){
    if(r.opt[i]>0){cover.num1 = c(cover.num1,length(which(dist.matM[i,]<r.opt[i])))}
    if(r.opt[i]==0 & min(dist.matM[i,])==0){cover.num1 = c(cover.num1,0.5)}
    if(r.opt[i]==0 & min(dist.matM[i,])>0){cover.num1 = c(cover.num1,0)}
  }
  
  MDS.index = which(cover.num1==max(cover.num1))   # The objects with most covering numbers are selected for an element of MDS at current step
  if(length(MDS.index)>1){
    MDS.index = MDS.index[which.max(ks.seq[MDS.index])[1]]
  } # Break ties with K-S statistic
    
  MDS = rbind(MDS,c(data[MDS.index,],r.opt[MDS.index],cover.num[MDS.index],ks.seq[MDS.index]))
  cover.MDS = c(cover.MDS,cover[MDS.index])
  rm = which(dist.matM[MDS.index,]<r.opt[MDS.index]) #objects covered at the current step (need to be removed)
  if(r.opt[MDS.index]==0){
    rm = which(dist.matM[MDS.index,]==0)
  }
  #r.opt1 = r.opt1[-MDS.index]
  #ks.seq1 = ks.seq1[-MDS.index]
  if(length(rm)==dim(dist.matM)[2]){
    len0=0
  } else {
    dist.matM = as.matrix(dist.matM[,-rm])
    len0 = dim(dist.matM)[2]
  }
}
colnames(MDS) = c("X","Y","r","coverNum","K-S")

for(i in 1:length(MDS[,1])){
  draw.circle(MDS[i,1],MDS[i,2],MDS[i,3],nv=1000,lwd=0.01)
}# draw covering balls



##### Build intersection graph for the approx MDS #####
lenM = length(MDS[,1]) #    # of the approx MDS

# distance matrix of the appro MDS
dist.MDS = NULL
for(a in 1:lenM){
  dist = NULL
  for(b in 1:lenM){
    dist.temp = sqrt((MDS[a,1]-MDS[b,1])^2+(MDS[a,2]-MDS[b,2])^2) #dist.temp: the distance between point a and b
    dist = c(dist,dist.temp) # dist is a vector representing the distances of $a$ with other points in the MDS 
  }
  dist.MDS = rbind(dist.MDS,dist)  # determine the distance between the current point and all others.
}
colnames(dist.MDS) = c(1:lenM)
rownames(dist.MDS) = c(1:lenM)

inter.num.MDS = NULL # record the intersection numbers of MDS
for(a in 1:lenM){
  count = 0
  for(b in 1:lenM){
    if(length(intersect(cover.MDS[[a]],cover.MDS[[b]]))>0){count = count + 1}
  }
  inter.num.MDS = c(inter.num.MDS,count)
}
print(inter.num.MDS)

MDS = cbind(MDS, inter.num.MDS)


###### Find the appro MDS of the intersect graph using algo 3  #####
MDS.I = NULL # store the elements of the approximate MDS for the intersect graph
MDS.r = MDS
dist.MDSM = dist.MDS #the distance matrix used for finding the MDS of the intersect graph
cover.MDS.I = cover.MDS 
len = dim(MDS.r)[1]

while(len>0){
  MDS.I = rbind(MDS.I,MDS.r[1,])
  if(len>1){
    rm = 1
    for(i in 2:len){
      if(length(intersect(cover.MDS.I[[1]],cover.MDS.I[[i]]))>0){rm = c(rm,i)}
    }
  }
  MDS.r = as.matrix(MDS.r[-rm,])
  if(len == 2){MDS.r = t(MDS.r)}
  cover.MDS.I = cover.MDS.I[-rm]
  len = dim(MDS.r)[1]
}


colnames(MDS) = c("X","Y","r","coverNum","K-S")

plot(data[1:100,1],data[1:100,2],xlab = "X", ylab = "Y",pch=20,ylim = c(2,20),xlim=c(-3,22),col="red",main="The Covering Balls of An Approx MDS (alg 3) of the intersect graph")
points(data[101:220,1],data[101:220,2],xlab = "X", ylab = "Y",pch=20,col="darkgreen")
points(data[221:370,1],data[221:370,2],xlab = "X", ylab = "Y",pch=20,col="blue")
points(data[371:373,1],data[371:373,2],xlab = "X", ylab = "Y",pch=3)
text(17,3,expression(paste(delta, " = 5")),cex=1.5)
legend("topleft",                                  
       legend=c("Cluster1","Cluster2","Cluster3","Outliers"),        
       col=c("red","darkgreen","blue","black"),                 
       pch=c(20,20,20,3),cex=1.5)
for(i in 1:length(MDS.I[,1])){
  draw.circle(MDS.I[i,1],MDS.I[i,2],MDS.I[i,3],nv=1000,lwd=0.01)
}# draw covering balls


###### Calculate the silhouette index and find the optimal clusters.
MDS.S = MDS.I[which(MDS.I[,3] != 0),]
sil.ave = NULL #record the  average silhouette index for each newly added cluster
len = length(MDS.S[,1])

for(a in 2:len){
  #assign each object to a cluster
  CluLabel = rep(0,length(data[,1]))
  for(i in 1:length(data[,1])){
    dist.temp = NULL
    for(j in 1:a){
      dist.ttemp = sqrt((data[i,1]-MDS.S[j,1])^2 + (data[i,2]-MDS.S[j,2])^2)/MDS.S[j,3]
      dist.temp = c(dist.temp, dist.ttemp)
    }
    CluLabel[i] = which.min(dist.temp)[1]
  }
  
  #Calculate the silhouette index for each objects
  sil = NULL
  for(i in 1:length(data[,1])){
    
       #compute a(i)
    Label.temp = CluLabel[i] # the label of the current object
    index.temp = which(CluLabel == Label.temp)
    index.temp = index.temp[-which(index.temp==i)]
    a.i = mean(dist.mat[i,index.temp])
    if(is.na(a.i)){a.i = 0}
    
    #compute b(i)
    b.i = NULL
    for(j in 1:a){
      if(j==Label.temp){next}
      index.temp = which(CluLabel == j)
      b.i = c(b.i, mean(dist.mat[i,index.temp]))
    }
    b.i = min(b.i)
    sil.i = (b.i-a.i)/max(a.i,b.i)
    sil = c(sil, sil.i)
  }
  sil.ave = c(sil.ave, mean(sil,na.rm = T))
}

index.temp = which.max(sil.ave) + 1
ClusCenter = MDS.I[1:index.temp,] ## The center of the clusters (the end) ##





plot(data[1:100,1],data[1:100,2],xlab = "X", ylab = "Y",pch=20,ylim = c(2,20),xlim=c(-3,22),col="red",main="The Covering Balls of An Approx MDS")
points(data[101:220,1],data[101:220,2],xlab = "X", ylab = "Y",pch=20,col="darkgreen")
points(data[221:370,1],data[221:370,2],xlab = "X", ylab = "Y",pch=20,col="blue")
points(data[371:373,1],data[371:373,2],xlab = "X", ylab = "Y",pch=3)
text(17,3,expression(paste(delta, " = 5")),cex=1.5)
legend("topleft",                                  
       legend=c("Cluster1","Cluster2","Cluster3","Outliers"),        
       col=c("red","darkgreen","blue","black"),                 
       pch=c(20,20,20,3),cex=1.5)

for(i in 1:length(dom[,1])){
  if(i == 1){draw.circle(dom[i,1],dom[i,2],dom[i,3],nv=1000,lty=1)}
  else {draw.circle(dom[i,1],dom[i,2],dom[i,3],nv=1000,lty=2)}
}





#Try delta = 7
plot(data[1:100,1],data[1:100,2],xlab = "X", ylab = "Y",pch=20,ylim = c(2,20),xlim=c(-3,22),col="red")
points(data[101:220,1],data[101:220,2],xlab = "X", ylab = "Y",pch=20,col="darkgreen")
points(data[221:370,1],data[221:370,2],xlab = "X", ylab = "Y",pch=20,col="blue")
points(data[371:373,1],data[371:373,2],xlab = "X", ylab = "Y",pch=3)
text(17,3,expression(paste(delta, " = 7")),cex=1.5)
legend("topleft",                                  
       legend=c("Cluster1","Cluster2","Cluster3","Outliers"),        
       col=c("red","darkgreen","blue","black"),                 
       pch=c(20,20,20,3),cex=1.5)      
#Set density parameter to be 10
delta = 7
d = 2
dist.all = NULL
#determine the radius of the covering ball for each points
cover.num0 = NULL
r.opt = NULL

for(a in 1:length(data[,1])){
  
  dist = NULL
  for(b in 1:length(data[,1])){
    dist.temp = sqrt((data[a,1]-data[b,1])^2+(data[a,2]-data[b,2])^2)
    dist = c(dist,dist.temp)
  }
  dist.all = rbind(dist.all,dist)  # determine the distance between the current point and all others.
  
  r.seq = seq(0.01,10.01,by=0.01)
  ks.seq = NULL
  for(i in 1:length(r.seq)){
    r.temp = r.seq[i]
    ks.temp = length(which(dist<r.temp))-delta*(r.temp^d)
    ks.seq =c(ks.seq,ks.temp)
  }
  cover.num0 = c(cover.num0,length(which(dist<r.seq[which.max(ks.seq)])))
  r.opt = c(r.opt,r.seq[which.max(ks.seq)]) # calculate the KS-statistics for the radii from 0.1 to 10 and find the maximum one.
  draw.circle(data[a,1],data[a,2],r.seq[which.max(ks.seq)],nv=1000,lwd=0.01)
}






plot(data[1:100,1],data[1:100,2],xlab = "X", ylab = "Y",pch=20,ylim = c(2,20),xlim=c(-3,22),col="red",main="The Covering Balls of An Approx MDS")
points(data[101:220,1],data[101:220,2],xlab = "X", ylab = "Y",pch=20,col="darkgreen")
points(data[221:370,1],data[221:370,2],xlab = "X", ylab = "Y",pch=20,col="blue")
points(data[371:373,1],data[371:373,2],xlab = "X", ylab = "Y",pch=3)
text(17,3,expression(paste(delta, " = 7")),cex=1.5)
legend("topleft",                                  
       legend=c("Cluster1","Cluster2","Cluster3","Outliers"),        
       col=c("red","darkgreen","blue","black"),                 
       pch=c(20,20,20,3),cex=1.5)  

# Finding an approximate minimal dominating set
dom = NULL #store the elements in approximate minimal dominating set
data0 = as.matrix(data)
r.opt0 = r.opt
dist.all0 = dist.all
while(length(data0[,1])>0){
  #re-calculate the covering number of each balls
  cover.num = NULL
  for(i in 1:length(data0[,1])){
    cover.num = c(cover.num,length(which(dist.all0[i,]<r.opt0[i])))
  }
  
  dom.ele = which.max(cover.num)[1]
  dom = rbind(dom,c(data0[dom.ele,],r.opt0[dom.ele]))
  rm = which(dist.all0[dom.ele,]<r.opt0[dom.ele])
  data0 = matrix(data0[-rm,],ncol=2)
  r.opt0 = r.opt0[-rm]
  dist.all0 = as.matrix(dist.all0[-rm,-rm])
}

for(i in 1:length(dom[,1])){
  draw.circle(dom[i,1],dom[i,2],dom[i,3],nv=1000,lwd=1)
}



# Build intersection graph for the approx MDS
inter.num = NULL
for(a in 1:length(dom[,1])){
  count = 0
  for(b in 1:length(dom[,1])){
    dist = sqrt((dom[a,1]-dom[b,1])^2 + (dom[a,2]-dom[b,2])^2)
    if(dist < (dom[a,3]+dom[b,3])){count = count + 1}
  }
  inter.num = c(inter.num,count)
}

dom.index = NULL
for(i in 1:length(dom[,1])){
  index.temp = which(data[,1]==dom[i,1])
  dom.index = c(dom.index,index.temp)
}

dom = cbind(dom,inter.num,cover.num0[dom.index])
colnames(dom)=c("x","y","r","inter.num","cover.num")
dom[order(-dom[,4],-dom[,5]),] #Sorting



#Obtain approximate MDS for the intersection graph


#Will incorporate silhouse index for cluster selection in the following few days.

