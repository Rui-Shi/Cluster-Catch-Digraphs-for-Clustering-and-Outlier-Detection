tt1 = Sys.time()
library(parallel)
library(doParallel)
library(MASS)
library(igraph)
source("G:/code_working_folder/Algo_Compare/LOF/LOF.R")
source("G:/code_working_folder/Algo_Compare/Real Datasets/RealData_Collection.R")
source("G:/code_working_folder/general functions/count.R")

cores = detectCores()

LB = 11 # Lower bound for MinPts
UB = 50 # Upper bound for MinPts
Thresh = 1.5 # threshhold for outliers



#### hepatitis data set ####
df = hepatitis
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels

score = LOF(X, LB, UB)

results = count_scores2(Y, score, Thresh)

print(paste("LOF x hepatitis: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### lymphography data set ####
df = lymphography
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels

score = LOF(X, LB, UB)

results = count_scores2(Y, score, Thresh)

print(paste("LOF x lymphography: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### glass data set ####
df = glass
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels

score = LOF(X, LB, UB)

results = count_scores2(Y, score, Thresh)

print(paste("LOF x glass: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### WBC data set ####
df = WBC
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels

score = LOF(X, LB, UB)

results = count_scores2(Y, score, Thresh)

print(paste("LOF x WBC: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### vertebral data set ####
#Thresh = 1.1 # threshhold for outliers
df = vertebral
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels

score = LOF(X, LB, UB)

results = count_scores2(Y, score, 1.1)

print(paste("LOF x vertebral: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### ecoli  data set ####
df = ecoli
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels

score = LOF(X, LB, UB)

results = count_scores2(Y, score, Thresh)

print(paste("LOF x ecoli: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### stamps  data set ####
df = stamps
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels

score = LOF(X, LB, UB)

results = count_scores2(Y, score, 1.2)

print(paste("LOF x stamps: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### WDBC  data set ####
df = WDBC
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels

score = LOF(X, LB, UB)

results = count_scores2(Y, score, Thresh)

print(paste("LOF x WDBC: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### pima data set ####
df = pima
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels

score = LOF(X, LB, UB)

results = count_scores2(Y, score, 1.2)

print(paste("LOF x pima: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### shuffle data set ####
df = shuffle
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels

score = LOF(X, LB, UB)

results = count_scores2(Y, score, Thresh)

print(paste("LOF x shuffle: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### vowels data set ####
#Thresh = 1.2 # threshhold for outliers
df = vowels
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels

score = LOF(X, LB, UB)

results = count_scores2(Y, score, 1.2)

print(paste("LOF x vowels: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### PenDigits  data set ####
# Thresh = 1.25 # threshhold for outliers
df = PenDigits
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels

score = LOF(X, LB, UB)

results = count_scores2(Y, score, 1.2)

print(paste("LOF x PenDigits: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### waveform  data set ####
#Thresh = 1.2 # threshhold for outliers
df = waveform
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels

score = LOF(X, LB, UB)

results = count_scores2(Y, score, 1.2)

print(paste("LOF x waveform: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### thyroid  data set ####
#Thresh = 1.2 # threshhold for outliers
df = thyroid
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels

score = LOF(X, LB, UB)

results = count_scores2(Y, score, 1.2)

print(paste("LOF x thyroid: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### wilt  data set ####
#Thresh = 1.2 # threshhold for outliers
df = wilt
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels

score = LOF(X, LB, UB)

results = count_scores2(Y, score, 1.2)

print(paste("LOF x wilt: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



#### pageblocks data set ####
df = pageblocks
n = dim(df)[1]
d = dim(df)[2]-1
X = df[,1:d] # predictors
Y = df[,d+1] # ground-true labels

score = LOF(X, LB, UB)

results = count_scores2(Y, score, 1.2)

print(paste("LOF x pageblocks: TPR is", results[1],",","TNR is", results[2],",","BA is", results[3],",","F2 is", results[4]))



tt2 = Sys.time()
tt2 - tt1