setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load library ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(missForest)
library(MSnbase)
library(pcaMethods)
library(VIM)
library(tictoc)  ## measure running time

# Load data ---------------------------------------------------------------

ori <- read.table("ni.3693-S2.txt",header = T,stringsAsFactors = F, quote = "",sep = "\t")

ori1 <- ori[,grep("LFQ.intensity_",colnames(ori))]

ori1 <- ori1[,grep("B\\.naive.*steady\\.state|T4\\.naive.*steady\\.state|T8\\.naive.*steady\\.state|MO.classical.*steady\\.state|B\\.naive.*activated|T4\\.naive.*activated|T8\\.naive.*activated|MO.classical.*activated",colnames(ori1))]

rownames(ori1) <- ori$Majority.protein.IDs
ori_anno <- ori[,1:5]
rownames(ori_anno) <- ori_anno$Majority.protein.IDs

ori2 <- ori1[rowSums(ori1 == 0) == 0,]  ## 3963 proteins
ori3 <- ori1[rowSums(ori1 == 0) < 31,]  ## 10081 proteins
sum(ori3 == 0)/10081/31  ## 27.07%

ori2 <- ori2[,c(16:31,1:15)]
ori2 <- log2(ori2)

ori3 <- ori3[,c(16:31,1:15)]
ori3[ori3 == 0] <- NA
ori3 <- log2(ori3)

source("functions_cell lineage.R")

# generate datasets with missing values -----------------------------------

## use MV = 0.1, 0.2, 0.3 and MNAR = 0.2, 0.5 0.8

mv <- c(rep(0.1,30),rep(0.2,30),rep(0.3,30))
mnar <- c(rep(c(rep(0.2,10),rep(0.5,10),rep(0.8,10)),3))
idx <- c(rep(1:10,9))

set.seed(88)
pick.seeds <- sample(1:10000,90)

cell_miss <- list()
for(i in 1:90){
  
  tmp <- addMiss(ori2, MV.rate = mv[i],MNAR.ratio = mnar[i], ini.seed = pick.seeds[i])
  tmp <- tmp[rowSums(is.na(tmp)) < ncol(tmp),]
  cell_miss[[i]] <- tmp
  
}

names(cell_miss) <- paste0("MV_",mv,"_MNAR_",mnar,"_",idx)

# run imputation and record run time --------------------------------------

time.table <- data.frame(matrix(NA,nrow = 90, ncol = 7))
rownames(time.table) <- names(cell_miss)
colnames(time.table) <- c("LOD","ND","kNN","LLS","RF","SVD","BPCA")

## LOD

cell_lod <- list()

for(i in 1:90){
  
  tictoc::tic("LOD")
  tmp <- cell_miss[[i]]
  tmp[is.na(tmp)] <- min(tmp,na.rm = T)
  tmp2 <- toc()
  time.table$LOD[i] <- as.numeric(tmp2$toc - tmp2$tic)
  cell_lod[[i]] <- tmp
}

## ND

cell_nd <- list()

for(i in 1:90){
  
  tictoc::tic("ND")
  cell_nd[[i]] <- normImp(cell_miss[[i]],width = 0.3,group = factor(c(rep("A",4),rep("B",4),rep("C",4),rep("D",4),rep("E",4),rep("F",4),rep("G",4),rep("H",3)),levels = c("A","B","C","D","E","F","G","H")),down.shift = 2.2,ori.seed = 666)
  tmp2 <- toc()
  time.table$ND[i] <- as.numeric(tmp2$toc - tmp2$tic)
  
}


## kNN, use k = 6

cell_knn <- list()

for(i in 1:90){
  
  tictoc::tic("kNN")
  tmp <- kNN(cell_miss[[i]],k = 6)
  tmp2 <- toc()
  time.table$kNN[i] <- as.numeric(tmp2$toc - tmp2$tic)
  
  rownames(tmp) <- rownames(cell_miss[[i]])
  cell_knn[[i]] <- tmp[,1:31]
  
  cat(paste("kNN completed",i,"datasets\n",collapse = " "))
}


## LLS

cell_lls <- list()

for(i in 1:90){
  
  tictoc::tic("LLS")
  tmp <- llsImpute(t(cell_miss[[i]]),allVariables = T, k = 150)
  tmp2 <- toc()
  time.table$LLS[i] <- as.numeric(tmp2$toc - tmp2$tic)
  
  cell_lls[[i]] <- as.data.frame(t(tmp@completeObs))
  
  cat(paste("LLS completed",i,"datasets\n",collapse = " "))
}


## RF

cell_rf <- list()

for(i in 1:90){
  
  tictoc::tic("RF")
  tmp <- missForest(cell_miss[[i]])
  tmp2 <- toc()
  time.table$RF[i] <- as.numeric(tmp2$toc - tmp2$tic)
  
  cell_rf[[i]] <- tmp$ximp
  
  cat(paste("RF completed",i,"datasets\n",collapse = " "))
}


## SVD, need to determine optimum nPCs!!

cell_svd <- list()

for(i in 1:90){
  
  tictoc::tic("SVD")
  tmp <- pca(cell_miss[[i]], method="svdImpute", nPcs=2, center = TRUE)
  tmp2 <- toc()
  time.table$SVD[i] <- as.numeric(tmp2$toc - tmp2$tic)
  
  cell_svd[[i]] <- as.data.frame(completeObs(tmp))
  
  cat(paste("SVD completed",i,"datasets\n",collapse = " "))
}


## BPCA, need to determine optimum nPCs!!

cell_bpca <- list()

for(i in 1:90){
  
  tictoc::tic("BPCA")
  tmp <- pca(cell_miss[[i]], method="bpca", nPcs=2)
  tmp2 <- toc()
  time.table$BPCA[i] <- as.numeric(tmp2$toc - tmp2$tic)
  
  cell_bpca[[i]] <- as.data.frame(completeObs(tmp))
  
  cat(paste("BPCA completed",i,"datasets\n",collapse = " "))
}



# save datasets -----------------------------------------------------------

save(ori2,ori3,ori_anno,cell_miss,cell_lod,cell_nd,cell_knn,cell_bpca,cell_svd,cell_rf,cell_lls,time.table, file = "201912_imputation data files_cell.RData")











