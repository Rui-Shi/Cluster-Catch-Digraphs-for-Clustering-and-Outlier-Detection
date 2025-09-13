tt1 = Sys.time()
library(parallel)
library(doParallel)
library(MASS)
library(igraph)
library(dbscan)
source("G:/code_working_folder/Algo_Compare/Isolation Forest/ISO.R")
source("G:/code_working_folder/Algo_Compare/Real Datasets/RealData_Collection.R")

ntrees = 1000; sample_size = 256; thresh = 0.55



#### hepatitis data set ####
df = hepatitis
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels
n0 = sum(Y==0) # number of outliers
n1 = n-n0

labels = ISO(X, ntrees, sample_size, thresh)

results = count_ISO2(n1, n0, labels)

print(paste("ISO x hepatitis: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### lymphography data set ####
df = lymphography
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels
n0 = sum(Y==0) # number of outliers
n1 = n-n0

labels = ISO(X, ntrees, sample_size, thresh)

results = count_ISO2(n1, n0, labels)

print(paste("ISO x lymphography: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### glass data set ####
df = glass
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels
n0 = sum(Y==0) # number of outliers
n1 = n-n0

labels = ISO(X, ntrees, sample_size, thresh)

results = count_ISO2(n1, n0, labels)

print(paste("ISO x glass: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### WBC data set ####
df = WBC
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels
n0 = sum(Y==0) # number of outliers
n1 = n-n0

labels = ISO(X, ntrees, sample_size, thresh)

results = count_ISO2(n1, n0, labels)

print(paste("ISO x WBC: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### vertebral data set ####
df = vertebral
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels
n0 = sum(Y==0) # number of outliers
n1 = n-n0

labels = ISO(X, ntrees, sample_size, thresh)

results = count_ISO2(n1, n0, labels)

print(paste("ISO x vertebral: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### ecoli data set ####
df = ecoli
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels
n0 = sum(Y==0) # number of outliers
n1 = n-n0

labels = ISO(X, ntrees, sample_size, threshold = 0.55)

results = count_ISO2(n1, n0, labels)

print(paste("ISO x ecoli: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### stamps data set ####
df = stamps
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels
n0 = sum(Y==0) # number of outliers
n1 = n-n0

labels = ISO(X, ntrees, sample_size, threshold = 0.55)

results = count_ISO2(n1, n0, labels)

print(paste("ISO x stamps: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### WDBC data set ####
df = WDBC
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels
n0 = sum(Y==0) # number of outliers
n1 = n-n0

labels = ISO(X, ntrees, sample_size, threshold = 0.55)

results = count_ISO2(n1, n0, labels)

print(paste("ISO x WDBC: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### pima data set ####
df = pima
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels
n0 = sum(Y==0) # number of outliers
n1 = n-n0

labels = ISO(X, ntrees, sample_size, threshold = 0.55)

results = count_ISO2(n1, n0, labels)

print(paste("ISO x pima: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))


#### shuffle data set ####
df = shuffle
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels
n0 = sum(Y==0) # number of outliers
n1 = n-n0

labels = ISO(X, ntrees, sample_size, threshold = 0.55)

results = count_ISO2(n1, n0, labels)

print(paste("ISO x shuffle: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### vowels data set ####
df = vowels
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels
n0 = sum(Y==0) # number of outliers
n1 = n-n0

labels = ISO(X, ntrees, sample_size, threshold = 0.55)

results = count_ISO2(n1, n0, labels)

print(paste("ISO x vowels: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### PenDigits data set ####
df = PenDigits
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels
n0 = sum(Y==0) # number of outliers
n1 = n-n0

labels = ISO(X, ntrees, sample_size, threshold = 0.55)

results = count_ISO2(n1, n0, labels)

print(paste("ISO x PenDigits: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### waveform data set ####
df = waveform
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels
n0 = sum(Y==0) # number of outliers
n1 = n-n0

labels = ISO(X, ntrees, sample_size, threshold = 0.55)

results = count_ISO2(n1, n0, labels)

print(paste("ISO x waveform: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### thyroid data set ####
df = thyroid
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels
n0 = sum(Y==0) # number of outliers
n1 = n-n0

labels = ISO(X, ntrees, sample_size, threshold = 0.55)

results = count_ISO2(n1, n0, labels)

print(paste("ISO x thyroid: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### wilt data set ####
df = wilt
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels
n0 = sum(Y==0) # number of outliers
n1 = n-n0

labels = ISO(X, ntrees, sample_size, threshold = 0.55)

results = count_ISO2(n1, n0, labels)

print(paste("ISO x wilt: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### pageblocks data set ####
df = pageblocks
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels
n0 = sum(Y==0) # number of outliers
n1 = n-n0

labels = ISO(X, ntrees, sample_size, threshold = 0.55)

results = count_ISO2(n1, n0, labels)

print(paste("ISO x pageblocks: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



tt2 = Sys.time()
tt2 - tt1