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
library(graphics)
library(gridExtra)
library(Hmisc)
library(FactoMineR)
library(DataExplorer)
library(factoextra)
library(expm)
library(adegraphics)
library(fpc)
library(cluster)

DataExplorer#Shift + cmd + c
setwd("/Users/JaviFerrando/Desktop/MVA-Project")

heart_disease = read.csv("data/heart.csv")
columns <- colnames(heart_disease)
columns[1] <- "age"
colnames(heart_disease) <- columns

# Find missing variables
which(is.na(heart_disease))

head(heart_disease)
describe(heart_disease)

#No missing values detected
#No need to perform any imputation


# Find missing variables

# if(which(is.na(heart_disease))){
#   heart_disease <- knnImputation(heart_disease, k = 1, scale = T)
# }

classVar <- lapply(heart_disease,class)   # class of each variable
factor_heart <- heart_disease
factor_heart$target <- as.factor(heart_disease$target)
factor_heart$sex <- as.factor(heart_disease$sex)
factor_heart$fbs <- as.factor(heart_disease$fbs)
factor_heart$exang <- as.factor(heart_disease$exang)
factor_heart$restecg <- as.factor(heart_disease$restecg)
factor_heart$thal <- as.factor(heart_disease$thal)
factor_heart$slope <- as.factor(heart_disease$slope)
factor_heart$cp <- as.factor(heart_disease$cp)
factor_heart$ca <- as.factor(heart_disease$ca)

#Outlier detection
############################
Moutlier(heart_disease[,-14], quantile = 0.975, plot = TRUE, tol=1e-36) #Doesn't work

pca_facto <- factor_heart[, sapply(factor_heart, class) != "factor"]
mout <- Moutlier(pca_facto, quantile = 0.975, plot = TRUE)
plot(mout$md,mout$rd,xlab='Classical Mahalanobis distance',ylab='Robust Mahalanobis distance')
abline(h = mout$cutoff, col="red")  # add cutoff line
abline(v = mout$cutoff, col="red")  # add cutoff line

#sort(round(mout$rd[mout$rd >= mout$cutoff],3),decreasing=T)

MRD_index_ordered <- order(mout$rd, decreasing=T)
MRD_index_ordered
round(mout$rd[MRD_index_ordered][1:5],3)

MRD_plot <- plot(mout$rd, 
                 pch="o", 
                 cex=1, 
                 main="Potential MRD outliers\n by Mahalanobis robust distance (MRD)",
                 ylab="MRD Rank")
#Local Outlier Factor
outlier.scores <- lofactor(heart_disease[,-14], k=5)
outlier.scores <- lofactor(pca_facto, k=5)
plot(density(outlier.scores),main='Distribution of individuals local outlier factor scores')

LOF_plot <- plot(outlier.scores, 
                 pch="o", 
                 cex=1, 
                 main="Potential LOF outliers\n by local outliers factor analysis (LOF-k=5)",
                 ylab="LOF Rank")
#LOF_plot_cutoff <- 0.5*(LOF_df[LOF_index_ordered[4],]$LOF_rank + LOF_df[LOF_index_ordered[5],]$LOF_rank)
abline(h = 1, col="red")  # add cutoff line

sort(outlier.scores, decreasing = TRUE)[1:4]

#Exploratory Data Analysis
#Density of heart presence/absence disease by age
g1 <- ggplot(data=heart_disease, aes(x=age, fill=as.factor(target)))+
  geom_density(alpha=.5)+
  ggtitle("Age") +
  scale_fill_manual(values = c('skyblue4', 'skyblue2'),name = "Disease", labels = c("Yes", "No"))

#Density of heart presence/absence disease by Max heart rate
g2 <- ggplot(data=heart_disease, aes(x=thalach, fill=as.factor(target)))+
  geom_density(alpha=.5)+
  ggtitle("Max Hear Rate") +
  scale_fill_manual(values = c('skyblue4', 'skyblue2'),name = "Disease", labels = c("Yes", "No"))

#Density of heart presence/absence disease by sex
g3 <- ggplot(data=heart_disease, aes(x=sex, fill=as.factor(target)))+
      geom_bar(alpha=.5, color="black")+
      ggtitle("Sex") +
      scale_fill_manual(values = c('skyblue4', 'skyblue2'),name = "Disease", labels = c("Yes", "No"))

#Density of heart presence/absence disease by chest type
g4 <- ggplot(data=heart_disease, aes(x=cp, fill=as.factor(target)))+
  geom_bar(alpha=.5, color="black")+
  ggtitle("Chest Pain type") +
  scale_fill_manual(values = c('skyblue4', 'skyblue2'),name = "Disease", labels = c("Yes", "No"))

grid.arrange(g1, g2, g3, g4, ncol = 2)

plot_correlation(heart_disease)

#PCA
pca_facto$sex <- factor_heart$sex
pca_facto$ca <- factor_heart$ca
pca_facto$disease <- heart_disease$target

pca_facto$disease[pca_facto$disease==0] <- "Yes"
pca_facto$disease[pca_facto$disease==1] <- "No"


pca_facto_heart <- PCA(pca_facto, quali.sup = 6, scale.unit = TRUE,  graph = TRUE)

#Screeplot
fviz_screeplot(pca_facto_heart, addlabels = FALSE)
eigen_values <- pca_facto_heart$eig[,1]
plot(eigen_values, type="o", main="Screeplot", 
     xlab='Dimension', ylab='Eigenvalue', col='blue')


#Represented in Rp
#quali.sup -> Every modality is the centroide of the respective individuals having chosen that modality
fviz_pca_ind(pca_facto_heart, habillage = 6, geom = "point", label="quali",addEllipses =TRUE, ellipse.level = 0.68)#co l.ind='cos2'
plot.PCA(pca_facto_heart, quali.sup = c(6:7), scale.unit = TRUE,choix = 'ind',label="quali")

#Represented in Rn
#Projection of variables, show correlation between principal components
fviz_pca_var(pca_facto_heart, geom = c("arrow", "text"), col.var = "cos2")#By quality of representation cos2


#Flipping
scaled.pca_facto <- scale(pca_facto)
n_row_X <- nrow(scaled.pca_facto)
matrix_N <- diag(x = 1/n_row_X, n_row_X, n_row_X, names = TRUE)
co_X <- t(scaled.pca_facto) %*% matrix_N %*% scaled.pca_facto# same as cor(scaled.pca_facto)

#R^p
U <- eigen(co_X)$vectors[,1:2] #U -> eigenvectors in R^p
proj_indiv <- scaled.pca_facto %*% U
plot(proj_indiv)

proj_indiv <- pca_facto_heart$ind$coord[,1:2]



#R^n
V <- eigen(sqrtm(matrix_N)%*%scaled.pca_facto%*%t(scaled.pca_facto)%*%sqrtm(matrix_N))$vectors #V -> eigenvectors in R^n
eigen(sqrtm(matrix_N)%*%scaled.pca_facto%*%t(scaled.pca_facto)%*%sqrtm(matrix_N))$values

proj_var <- t(scaled.pca_facto)%*%sqrtm(matrix_N)%*%V
plot(proj_var)
s.arrow(proj_var)


#Clustering
hc_ward = hclust(dist(proj_indiv),method = "ward.D")
plot(hc_ward, main= "HC using Ward Agglomeration method", xlab="",sub="",cex=.9)
abline(h=60)
rect.hclust(hc_ward, k = 3, border = 2:6)

#Association of individuals to clusters
classes <- cutree(hc_ward, h=50) #Depending on the height, number of clusters is chosen
plotcluster(proj_indiv, classes,main="Projections of individuals in Hirechical Clustering of 3 classes")

get_centroids <- function(classes, n_classes){
  centroids <- NULL
  for(k in 1:n_classes){
    centroids <- rbind(centroids, colMeans(proj_indiv[classes == k, , drop = FALSE]))
  }
  return(centroids)
}
centroids <- get_centroids(classes, 3)

k_mean <- kmeans(proj_indiv, centroids)
plotcluster(proj_indiv, k_mean$cluster,main="Projections of individuals in a K-means cluster of 3 classes")


cal_idx_before <- calinhara(proj_indiv,classes,cn=max(classes))
cal_idx_after <- calinhara(proj_indiv,k_mean$cluster,cn=max(k_mean$cluster))

print(cal_idx_before)
print(cal_idx_after)

#Calinski-Harabasz index
Calinski_Harabassza <- function (projections, hc, kind, n_classes){
  classes <- cutree(hc, k=n_classes)
  centroids <- get_centroids(classes, n_classes)
  if(kind=='hc'){
    index <- calinhara(proj_indiv,classes,cn=max(classes))
  }
  if(kind=='kmeans'){
    kmeans_classes <- kmeans(proj_indiv, centers = centroids)$cluster
    index <-calinhara(proj_indiv,kmeans_classes,cn=max(kmeans_classes)) 
  }
  return(index)
}
get_indexes <- function(until, kind){
  indexes <- c()
  for (n_classes in 2:until){
    indexes <- c(indexes, Calinski_Harabassza(proj_indiv, hc_ward, kind, n_classes))
  }  
  return(indexes)
}

####
indexes_before <- get_indexes(10, 'hc')
plot(indexes_before, type = "o", xlab = 'Number of classes', ylab = 'Calinsky index value'
     , main = 'Index before consolidation', col = 'firebrick1', xaxt
     = "n")
axis(1, at=1:9, labels = c(2, 3, 4, 5, 6, 7,8,9,10))

###
indexes_after <- get_indexes(10, 'kmeans')
plot(indexes_after, type = "o", xlab = 'Number of classes', ylab = 'Calinsky index value'
     , main = 'Index after consolidation', col = 'firebrick1', xaxt
     = "n")
axis(1, at=1:9, labels = c(2, 3, 4, 5, 6, 7,8,9,10))  

#Multiple Correspondence Analysis




