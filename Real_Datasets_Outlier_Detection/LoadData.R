setwd("G:/code_working_folder/Algo_Compare_OD/Real_Datasets")
library(R.matlab)
library(foreign)

glass = readMat("glass.mat")

vertebral = readMat("vertebral.mat")

vowels = readMat("vowels.mat")

Ecoli = readMat('Ecoli.mat')

PenDigits = read.arff('PenDigits_withoutdupl_norm_v01.arff')

Shuttle = read.arff('Shuttle_withoutdupl_norm_v01.arff')

Wilt = read.arff("Wilt_withoutdupl_norm_02_v01.arff")

Wave = read.arff('Waveform_withoutdupl_v01.arff')

seismic = read.arff('seismic-bumps.arff')

thyroid = readMat('thyroid.mat')
