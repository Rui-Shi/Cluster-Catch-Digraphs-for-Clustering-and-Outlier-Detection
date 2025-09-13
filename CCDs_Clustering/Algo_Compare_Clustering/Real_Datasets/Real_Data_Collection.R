library(R.matlab)
library(foreign)
library(dplyr)
library(readxl)
setwd("/media/rui/exNVME/code_working_folder/Algo_Compare_Clustering/Real_Datasets")
# setwd("G:/code_working_folder/Algo_Compare_Clustering/Real_Datasets")

# a robust version of normalization
scale_R = function(x){
  M = median(x)
  madn = mad(x)
  return((x-M)/madn)
}

##### iris #####
iris = read.csv("iris.data", header=F)
X = iris[, -ncol(iris)]

# normalization
X = apply(X, 2, scale)

# remove the observations with NA values
X = X[!rowSums(is.na(X)) > 0, ] 

Y = iris[,ncol(iris)]
Y = as.numeric(as.factor(Y))
iris = cbind(X, Y)

# remove duplicated rows
iris = as.matrix(distinct(as.data.frame(iris)))

# Check if there are NV values
anyNA(iris)



##### seeds #####
seeds = read.table("seeds_dataset.txt")

X = seeds[,-ncol(seeds)]

# normalization
X = apply(X, 2, scale)

# remove the observations with NA values
X = X[!rowSums(is.na(X)) > 0, ] 

Y = seeds[,ncol(seeds)]
seeds = cbind(X, Y)

# remove duplicated rows
seeds = as.matrix(distinct(as.data.frame(seeds)))

# Check if there are NV values
anyNA(seeds)



##### whole sale #####
wholesale =  read.csv("Wholesale customers data.csv")

X = wholesale[, -c(1,2)]
X = as.matrix(as.data.frame(lapply(X, as.numeric)))

# normalization
X = apply(X, 2, scale)

# remove the observations with NA values
X = X[!rowSums(is.na(X)) > 0, ] 

Y = wholesale[ ,2]
Y = as.numeric(as.factor(Y))
wholesale = cbind(X, Y)

# remove duplicated rows
wholesale = as.matrix(distinct(as.data.frame(wholesale)))

# Check if there are NA values
anyNA(wholesale)



##### breast cancer #####
breast_cancer = read.csv("breast-cancer-wisconsin.data", header=F)

breast_cancer = breast_cancer[rowSums(breast_cancer == "?") == 0, ] 

X = breast_cancer[, 2:(ncol(breast_cancer)-1)]
X = as.matrix(as.data.frame(lapply(X, as.numeric)))

# normalization
X = apply(X, 2, scale)

# remove the observations with NA values
X = X[!rowSums(is.na(X)) > 0, ] 

Y = breast_cancer[,ncol(breast_cancer)]
breast_cancer = cbind(X, Y)

# remove duplicated rows
breast_cancer = as.matrix(distinct(as.data.frame(breast_cancer)))

# Check if there are NV values
anyNA(breast_cancer)



##### aggregation #####
aggregation = read.table("aggregation.txt", header=FALSE)

X = aggregation[,-ncol(aggregation)]

# normalization
# X = apply(X, 2, scale)

# remove the observations with NA values
X = X[!rowSums(is.na(X)) > 0, ] 

Y = aggregation[,ncol(aggregation)]
aggregation = cbind(X, Y)

# remove duplicated rows
aggregation = as.matrix(distinct(as.data.frame(aggregation)))

# Check if there are NV values
anyNA(aggregation)



##### asymmetric #####
asymmetric = read.table("asymmetric.txt", header=FALSE)
Y = read.table("asymmetric_label.txt", header=FALSE)

X = asymmetric

# normalization
X = apply(X, 2, scale)

# remove the observations with NA values
X = X[!rowSums(is.na(X)) > 0, ] 

asymmetric = cbind(X, Y)

# remove duplicated rows
asymmetric = as.matrix(distinct(as.data.frame(asymmetric)))

# Check if there are NV values
anyNA(asymmetric)

# plot(asymmetric[,1], asymmetric[,2], pch=20, asp=1, xaxt='n', yaxt='n', ylab="", xlab="")


##### R15 #####
R15 = read.table("R15.txt", header=FALSE)

X = R15[,-ncol(R15)]

# normalization
# X = apply(X, 2, scale)

# remove the observations with NA values
X = X[!rowSums(is.na(X)) > 0, ] 

Y = R15[,ncol(R15)]
R15 = cbind(X, Y)

# remove duplicated rows
R15 = as.matrix(distinct(as.data.frame(R15)))

# Check if there are NV values
anyNA(R15)

# plot(R15[,1], R15[,2], pch=20, asp=1, xaxt='n', yaxt='n', ylab="", xlab="")



##### D31 #####
D31 = read.table("D31.txt", header=FALSE)

X = D31[,-ncol(D31)]

# normalization
# X = apply(X, 2, scale)

# remove the observations with NA values
X = X[!rowSums(is.na(X)) > 0, ] 

Y = D31[,ncol(D31)]
D31 = cbind(X, Y)

# remove duplicated rows
D31 = as.matrix(distinct(as.data.frame(D31)))

# Check if there are NV values
anyNA(D31)

# plot(D31[,1], D31[,2], pch=20, asp=1, xaxt='n', yaxt='n', ylab="", xlab="")


##### user_knowledge #####
user_knowledge = read_excel("user_knowlege.xls")

X = user_knowledge[,1:5]

# normalization
X = apply(X, 2, scale)

# remove the observations with NA values
X = X[!rowSums(is.na(X)) > 0, ]

# convert Y to numeric numbers
Y = user_knowledge$UNS
Y = as.numeric(as.factor(Y))

user_knowledge = cbind(X, Y)

# remove duplicated rows
user_knowledge = as.matrix(distinct(as.data.frame(user_knowledge)))

# Check if there are NA values
anyNA(user_knowledge)