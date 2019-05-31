# MVA PROJECT STEP (1)
#
#
# Pre-process of data. The student will perform a first summary of the
# data, and, eventually, detect errors, outliers and missing values and
# take the appropriate measures of correction. According to the
# problem and data, it may be necessary to perform a selection of
# variables (feature selection) and /or a derivation of new explanatory
# variables (feature extraction) if the problem needs it.
#

library(chemometrics)
library(DMwR)
library(mice)
library(missForest)
library(ggplot2)


#Shift + cmd + c
setwd("/Users/JaviFerrando/Desktop/MVA-Project")

heart_disease = read.csv("data/heart.csv")

# Find missing variables
which(is.na(heart_disease))

#No missing values detected
#No need to perform any imputation


# Find missing variables

# if(which(is.na(heart_disease))){
#   heart_disease <- knnImputation(heart_disease, k = 1, scale = T)
# }

classVar <- lapply(heart_disease,class)   # class of each variable

#Outlier detection
############################
Moutlier(heart_disease[,-14], quantile = 0.975, plot = TRUE, tol=1e-36) #Doesn't work


#Local Outlier Factor
outlier.scores <- lofactor(heart_disease[,-14], k=5)
plot(density(outlier.scores),main='Distribution of individuals local outlier factor scores')



