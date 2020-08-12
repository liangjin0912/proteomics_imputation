setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load library ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(missForest)
library(MSnbase)
library(pcaMethods)
library(VIM)
library(tictoc)  ## measure running time

# load data and functions -------------------------------------------------

ion <- read.xlsx("56 biological replicates_quantitative data.xlsx")
rownames(ion) <- ion$ID
ion <- ion[,c(-1,-2)]
ion <- ion[,1:32]

ion[ion <=  10] <- NA
ion <- ion[rowSums(is.na(ion)) == 0,]  ## 6689

source("functions_3 proteome.R")


# generate datasets with missing values -----------------------------------

## use MV = 0.1, 0.2, 0.3 and MNAR = 0.2, 0.5 0.8

mv <- c(rep(0.1,30),rep(0.2,30),rep(0.3,30))
mnar <- c(rep(c(rep(0.2,10),rep(0.5,10),rep(0.8,10)),3))
idx <- c(rep(1:10,9))

set.seed(233)
pick.seeds <- sample(1:10000,90)

ion_miss <- list()
for(i in 1:90){
  
  tmp <- addMiss(ion, MV.rate = mv[i],MNAR.ratio = mnar[i], ini.seed = pick.seeds[i])
  tmp <- tmp[rowSums(is.na(tmp)) < ncol(tmp),]
  ion_miss[[i]] <- tmp
  
}

names(ion_miss) <- paste0("MV_",mv,"_MNAR_",mnar,"_",idx)


# run imputation and record run time --------------------------------------

time.table <- data.frame(matrix(NA,nrow = 90, ncol = 7))
rownames(time.table) <- names(ion_miss)
colnames(time.table) <- c("LOD","ND","kNN","LLS","RF","SVD","BPCA")

## LOD

ion_lod <- list()

for(i in 1:90){
  
  tictoc::tic("LOD")
  tmp <- ion_miss[[i]]
  tmp[is.na(tmp)] <- min(tmp,na.rm = T)
  tmp2 <- toc()
  time.table$LOD[i] <- as.numeric(tmp2$toc - tmp2$tic)
  ion_lod[[i]] <- tmp
}

## ND

ion_nd <- list()

for(i in 1:90){
  
  tictoc::tic("ND")
  ion_nd[[i]] <- normImp(ion_miss[[i]],width = 0.3,group = factor(c(rep("A",8),rep("B",8),rep("C",8),rep("D",8)),levels = c("A","B","C","D")),down.shift = 2.2,ori.seed = 666)
  tmp2 <- toc()
  time.table$ND[i] <- as.numeric(tmp2$toc - tmp2$tic)
  
}


## kNN, use k = 6

ion_knn <- list()

for(i in 1:90){
  
  tictoc::tic("kNN")
  tmp <- kNN(ion_miss[[i]],k = 6)
  tmp2 <- toc()
  time.table$kNN[i] <- as.numeric(tmp2$toc - tmp2$tic)
  
  rownames(tmp) <- rownames(ion_miss[[i]])
  ion_knn[[i]] <- tmp[,1:32]
  
  cat(paste("kNN completed",i,"datasets\n",collapse = " "))
}


## LLS

ion_lls <- list()

for(i in 1:90){
  
  tictoc::tic("LLS")
  tmp <- llsImpute(t(ion_miss[[i]]),allVariables = T, k = 150)
  tmp2 <- toc()
  time.table$LLS[i] <- as.numeric(tmp2$toc - tmp2$tic)
  
  ion_lls[[i]] <- as.data.frame(t(tmp@completeObs))
  
  cat(paste("LLS completed",i,"datasets\n",collapse = " "))
}


## RF

ion_rf <- list()

for(i in 1:90){
  
  tictoc::tic("RF")
  tmp <- missForest(ion_miss[[i]])
  tmp2 <- toc()
  time.table$RF[i] <- as.numeric(tmp2$toc - tmp2$tic)
  
  ion_rf[[i]] <- tmp$ximp
  
  cat(paste("RF completed",i,"datasets\n",collapse = " "))
}


## SVD, need to determine optimum nPCs!!

ion_svd <- list()

for(i in 1:90){
  
  tictoc::tic("SVD")
  tmp <- pca(ion_miss[[i]], method="svdImpute", nPcs=2, center = TRUE)
  tmp2 <- toc()
  time.table$SVD[i] <- as.numeric(tmp2$toc - tmp2$tic)
  
  ion_svd[[i]] <- as.data.frame(completeObs(tmp))
  
  cat(paste("SVD completed",i,"datasets\n",collapse = " "))
}


## BPCA, need to determine optimum nPCs!!

ion_bpca <- list()

for(i in 1:90){
  
  tictoc::tic("BPCA")
  tmp <- pca(ion_miss[[i]], method="bpca", nPcs=2)
  tmp2 <- toc()
  time.table$BPCA[i] <- as.numeric(tmp2$toc - tmp2$tic)
  
  ion_bpca[[i]] <- as.data.frame(completeObs(tmp))
  
  cat(paste("BPCA completed",i,"datasets\n",collapse = " "))
}



# save datasets -----------------------------------------------------------

#save(ion,ion_miss,ion_lod,ion_nd,ion_knn,ion_bpca,ion_svd,ion_rf,ion_lls,time.table, file = "201911_imputation data files.RData")





