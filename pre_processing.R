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

#Shift + cmd + c
setwd("/Users/JaviFerrando/Desktop/MVA-Project")

heart_disease = read.csv("data/heart.csv")


# Find missing variables

if(which(is.na(heart_disease))){
  heart_disease <- knnImputation(heart_disease, k = 1, scale = T)
}
