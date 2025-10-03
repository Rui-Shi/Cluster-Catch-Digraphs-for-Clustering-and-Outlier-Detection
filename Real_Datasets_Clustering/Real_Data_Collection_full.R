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

##### wine-red #####
wine_red = read.csv("winequality-red.csv", header=T, sep = ";")

X = wine_red[,-ncol(wine_red)]

# normalization
X = apply(X, 2, scale)

# remove the observations with NA values
X = X[!rowSums(is.na(X)) > 0, ] 

Y = wine_red[,ncol(wine_red)]
wine_red = cbind(X, Y)

# remove duplicated rows
wine_red = as.matrix(distinct(as.data.frame(wine_red)))

# Check if there are NV values
anyNA(wine_red)



##### wine-white #####
wine_white = read.csv("winequality-white.csv", header=T, sep = ";")

X = wine_white[,-ncol(wine_white)]

# normalization
X = apply(X, 2, scale_R)

# remove the observations with NA values
X = X[!rowSums(is.na(X)) > 0, ] 

Y = wine_white[,ncol(wine_white)]
wine_white = cbind(X, Y)

# remove duplicated rows
wine_white = as.matrix(distinct(as.data.frame(wine_white)))

# Check if there are NV values
anyNA(wine_white)



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



##### yeast #####
yeast = read.table("yeast.data", quote="\"", comment.char="")

X = yeast[, 2:(ncol(yeast)-1)]

# normalization
X = apply(X, 2, scale)

# remove the observations with NA values
X = X[!rowSums(is.na(X)) > 0, ] 

Y = yeast[,ncol(yeast)]
Y = as.numeric(as.factor(Y))
yeast = cbind(X, Y)

# remove duplicated rows
yeast = as.matrix(distinct(as.data.frame(yeast)))

# Check if there are NV values
anyNA(yeast)



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



##### transfusion #####
transfusion = read.csv("transfusion.data")

X = transfusion[, -ncol(transfusion)]

# normalization
X = apply(X, 2, scale)

# remove the observations with NA values
X = X[!rowSums(is.na(X)) > 0, ] 

Y = transfusion[,ncol(transfusion)]
Y = as.numeric(as.factor(Y))
transfusion = cbind(X, Y)

# remove duplicated rows
transfusion = as.matrix(distinct(as.data.frame(transfusion)))

# Check if there are NV values
anyNA(transfusion)



##### WebPhishing #####
phishing = read.arff('PhishingData.arff')

X = phishing[, -ncol(phishing)]
X = as.matrix(as.data.frame(lapply(X, as.numeric)))

# normalization
X = apply(X, 2, scale)

# remove the observations with NA values
X = X[!rowSums(is.na(X)) > 0, ] 

Y = phishing[,ncol(phishing)]
Y = as.numeric(as.factor(Y))
phishing = cbind(X, Y)

# remove duplicated rows
phishing = as.matrix(distinct(as.data.frame(phishing)))

# Check if there are NV values
anyNA(phishing)



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



##### glass #####
glass = read.table("glass.txt", sep = ",")

X = glass[,2:10]

# normalization
# X = apply(X, 2, scale)

# remove the observations with NA values
X = X[!rowSums(is.na(X)) > 0, ] 

Y = glass[,ncol(glass)]
glass = as.matrix(cbind(X, Y))

# remove duplicated rows
glass = as.matrix(distinct(as.data.frame(glass)))

# Check if there are NA values
anyNA(glass)



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



##### hepatitis #####
hepatitis = read.csv("hepatitis.data")
hepatitis = hepatitis[rowSums(hepatitis== "?") == 0, ]

X = hepatitis[, -ncol(hepatitis)]
X = as.matrix(as.data.frame(lapply(X, as.numeric)))

# normalization
X = apply(X, 2, scale)

# remove the observations with NA values
X = X[!rowSums(is.na(X)) > 0, ] 

Y = hepatitis[,ncol(hepatitis)]
Y = as.numeric(as.factor(Y))
hepatitis = cbind(X, Y)

# remove duplicated rows
hepatitis = as.matrix(distinct(as.data.frame(hepatitis)))

# Check if there are NA values
anyNA(hepatitis)



##### lymphography #####
lymphography = read.csv("lymphography.data")
# lymphography = lymphography[rowSums(lymphography== "?") == 0, ]

X = lymphography[, -ncol(lymphography)]
# X = as.matrix(as.data.frame(lapply(X, as.numeric)))

# normalization
X = apply(X, 2, scale)

# remove the observations with NA values
X = X[!rowSums(is.na(X)) > 0, ] 

Y = lymphography[,ncol(lymphography)]
Y = as.numeric(as.factor(Y))
lymphography = cbind(X, Y)

# remove duplicated rows
lymphography = as.matrix(distinct(as.data.frame(lymphography)))

# Check if there are NA values
anyNA(lymphography)



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