setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# load data and library ---------------------------------------------------

load("201911_imputation data files.RData")
source("functions_3 proteome.R")

library(dplyr)
library(missForest)
library(purrr)
library(ggplot2)

# generate summary tables -------------------------------------------------

sum_ion <- SumTable(ion)

sum_lod <- lapply(ion_lod,SumTable)
names(sum_lod) <- names(ion_lod)

sum_nd <- lapply(ion_nd,SumTable)
names(sum_nd) <- names(ion_nd)

sum_knn <- lapply(ion_knn,SumTable)
names(sum_knn) <- names(ion_knn)

sum_lls <- lapply(ion_lls,SumTable)
names(sum_lls) <- names(ion_lls)

sum_rf <- lapply(ion_rf,SumTable)
names(sum_rf) <- names(ion_rf)

sum_svd <- lapply(ion_svd,SumTable)
names(sum_svd) <- names(ion_svd)

sum_bpca <- lapply(ion_bpca,SumTable)
names(sum_bpca) <- names(ion_bpca)

# generate result tables --------------------------------------------------

## lod

res_lod <- data.frame(matrix(NA,nrow = 90, ncol = 17))
colnames(res_lod) <- c("time","NRMSE","MAE.B.A","MAE.C.A","MAE.D.A","NRMSE.B.A","NRMSE.C.A","NRMSE.D.A","TP.B.A","TP.C.A","TP.D.A","FP.B.A","FP.C.A","FP.D.A","FADR.B.A","FADR.C.A","FADR.D.A")
rownames(res_lod) <- rownames(time.table)
res_lod$time <- time.table$LOD

for(i in 1:90){
  
  res_lod$NRMSE[i] <- missForest::nrmse(ion_lod[[i]],ion_miss[[i]],ion[rownames(ion_miss[[i]]),])
  res_lod$MAE.B.A[i] <- mean(abs(sum_lod[[i]]$FC.B.A - sum_lod[[i]]$ori.FC.B.A))
  res_lod$MAE.C.A[i] <- mean(abs(sum_lod[[i]]$FC.C.A - sum_lod[[i]]$ori.FC.C.A))
  res_lod$MAE.D.A[i] <- mean(abs(sum_lod[[i]]$FC.D.A - sum_lod[[i]]$ori.FC.D.A))
  res_lod$NRMSE.B.A[i] <- cal.nrmse(sum_lod[[i]]$ori.FC.B.A,sum_lod[[i]]$FC.B.A)
  res_lod$NRMSE.C.A[i] <- cal.nrmse(sum_lod[[i]]$ori.FC.C.A,sum_lod[[i]]$FC.C.A)
  res_lod$NRMSE.D.A[i] <- cal.nrmse(sum_lod[[i]]$ori.FC.D.A,sum_lod[[i]]$FC.D.A)
  res_lod$TP.B.A[i] <- sum(sum_lod[[i]]$adj.p.B.A < 0.05 & sum_lod[[i]]$species == "ECOLI")
  res_lod$TP.C.A[i] <- sum(sum_lod[[i]]$adj.p.C.A < 0.05 & sum_lod[[i]]$species == "ECOLI")
  res_lod$TP.D.A[i] <- sum(sum_lod[[i]]$adj.p.D.A < 0.05 & sum_lod[[i]]$species == "ECOLI")
  res_lod$FP.B.A[i] <- sum(sum_lod[[i]]$adj.p.B.A < 0.05 & sum_lod[[i]]$species == "HUMAN")
  res_lod$FP.C.A[i] <- sum(sum_lod[[i]]$adj.p.C.A < 0.05 & sum_lod[[i]]$species == "HUMAN")
  res_lod$FP.D.A[i] <- sum(sum_lod[[i]]$adj.p.D.A < 0.05 & sum_lod[[i]]$species == "HUMAN")

}

res_lod$FADR.B.A <- res_lod$FP.B.A/(res_lod$TP.B.A+res_lod$FP.B.A)
res_lod$FADR.C.A <- res_lod$FP.C.A/(res_lod$TP.C.A+res_lod$FP.C.A)
res_lod$FADR.D.A <- res_lod$FP.D.A/(res_lod$TP.D.A+res_lod$FP.D.A)

## nd

res_nd <- data.frame(matrix(NA,nrow = 90, ncol = 17))
colnames(res_nd) <- c("time","NRMSE","MAE.B.A","MAE.C.A","MAE.D.A","NRMSE.B.A","NRMSE.C.A","NRMSE.D.A","TP.B.A","TP.C.A","TP.D.A","FP.B.A","FP.C.A","FP.D.A","FADR.B.A","FADR.C.A","FADR.D.A")
rownames(res_nd) <- rownames(time.table)
res_nd$time <- time.table$ND

for(i in 1:90){
  
  res_nd$NRMSE[i] <- missForest::nrmse(ion_nd[[i]],ion_miss[[i]],ion[rownames(ion_miss[[i]]),])
  res_nd$MAE.B.A[i] <- mean(abs(sum_nd[[i]]$FC.B.A - sum_nd[[i]]$ori.FC.B.A))
  res_nd$MAE.C.A[i] <- mean(abs(sum_nd[[i]]$FC.C.A - sum_nd[[i]]$ori.FC.C.A))
  res_nd$MAE.D.A[i] <- mean(abs(sum_nd[[i]]$FC.D.A - sum_nd[[i]]$ori.FC.D.A))
  res_nd$NRMSE.B.A[i] <- cal.nrmse(sum_nd[[i]]$ori.FC.B.A,sum_nd[[i]]$FC.B.A)
  res_nd$NRMSE.C.A[i] <- cal.nrmse(sum_nd[[i]]$ori.FC.C.A,sum_nd[[i]]$FC.C.A)
  res_nd$NRMSE.D.A[i] <- cal.nrmse(sum_nd[[i]]$ori.FC.D.A,sum_nd[[i]]$FC.D.A)
  res_nd$TP.B.A[i] <- sum(sum_nd[[i]]$adj.p.B.A < 0.05 & sum_nd[[i]]$species == "ECOLI")
  res_nd$TP.C.A[i] <- sum(sum_nd[[i]]$adj.p.C.A < 0.05 & sum_nd[[i]]$species == "ECOLI")
  res_nd$TP.D.A[i] <- sum(sum_nd[[i]]$adj.p.D.A < 0.05 & sum_nd[[i]]$species == "ECOLI")
  res_nd$FP.B.A[i] <- sum(sum_nd[[i]]$adj.p.B.A < 0.05 & sum_nd[[i]]$species == "HUMAN")
  res_nd$FP.C.A[i] <- sum(sum_nd[[i]]$adj.p.C.A < 0.05 & sum_nd[[i]]$species == "HUMAN")
  res_nd$FP.D.A[i] <- sum(sum_nd[[i]]$adj.p.D.A < 0.05 & sum_nd[[i]]$species == "HUMAN")
  
}

res_nd$FADR.B.A <- res_nd$FP.B.A/(res_nd$TP.B.A+res_nd$FP.B.A)
res_nd$FADR.C.A <- res_nd$FP.C.A/(res_nd$TP.C.A+res_nd$FP.C.A)
res_nd$FADR.D.A <- res_nd$FP.D.A/(res_nd$TP.D.A+res_nd$FP.D.A)

## knn

res_knn <- data.frame(matrix(NA,nrow = 90, ncol = 17))
colnames(res_knn) <- c("time","NRMSE","MAE.B.A","MAE.C.A","MAE.D.A","NRMSE.B.A","NRMSE.C.A","NRMSE.D.A","TP.B.A","TP.C.A","TP.D.A","FP.B.A","FP.C.A","FP.D.A","FADR.B.A","FADR.C.A","FADR.D.A")
rownames(res_knn) <- rownames(time.table)
res_knn$time <- time.table$kNN

for(i in 1:90){
  
  res_knn$NRMSE[i] <- nrmse(ion_knn[[i]],ion_miss[[i]],ion[rownames(ion_miss[[i]]),])
  res_knn$MAE.B.A[i] <- mean(abs(sum_knn[[i]]$FC.B.A - sum_knn[[i]]$ori.FC.B.A))
  res_knn$MAE.C.A[i] <- mean(abs(sum_knn[[i]]$FC.C.A - sum_knn[[i]]$ori.FC.C.A))
  res_knn$MAE.D.A[i] <- mean(abs(sum_knn[[i]]$FC.D.A - sum_knn[[i]]$ori.FC.D.A))
  res_knn$NRMSE.B.A[i] <- cal.nrmse(sum_knn[[i]]$ori.FC.B.A,sum_knn[[i]]$FC.B.A)
  res_knn$NRMSE.C.A[i] <- cal.nrmse(sum_knn[[i]]$ori.FC.C.A,sum_knn[[i]]$FC.C.A)
  res_knn$NRMSE.D.A[i] <- cal.nrmse(sum_knn[[i]]$ori.FC.D.A,sum_knn[[i]]$FC.D.A)
  res_knn$TP.B.A[i] <- sum(sum_knn[[i]]$adj.p.B.A < 0.05 & sum_knn[[i]]$species == "ECOLI")
  res_knn$TP.C.A[i] <- sum(sum_knn[[i]]$adj.p.C.A < 0.05 & sum_knn[[i]]$species == "ECOLI")
  res_knn$TP.D.A[i] <- sum(sum_knn[[i]]$adj.p.D.A < 0.05 & sum_knn[[i]]$species == "ECOLI")
  res_knn$FP.B.A[i] <- sum(sum_knn[[i]]$adj.p.B.A < 0.05 & sum_knn[[i]]$species == "HUMAN")
  res_knn$FP.C.A[i] <- sum(sum_knn[[i]]$adj.p.C.A < 0.05 & sum_knn[[i]]$species == "HUMAN")
  res_knn$FP.D.A[i] <- sum(sum_knn[[i]]$adj.p.D.A < 0.05 & sum_knn[[i]]$species == "HUMAN")
  
}

res_knn$FADR.B.A <- res_knn$FP.B.A/(res_knn$TP.B.A+res_knn$FP.B.A)
res_knn$FADR.C.A <- res_knn$FP.C.A/(res_knn$TP.C.A+res_knn$FP.C.A)
res_knn$FADR.D.A <- res_knn$FP.D.A/(res_knn$TP.D.A+res_knn$FP.D.A)

## lls

res_lls <- data.frame(matrix(NA,nrow = 90, ncol = 17))
colnames(res_lls) <- c("time","NRMSE","MAE.B.A","MAE.C.A","MAE.D.A","NRMSE.B.A","NRMSE.C.A","NRMSE.D.A","TP.B.A","TP.C.A","TP.D.A","FP.B.A","FP.C.A","FP.D.A","FADR.B.A","FADR.C.A","FADR.D.A")
rownames(res_lls) <- rownames(time.table)
res_lls$time <- time.table$LLS

for(i in 1:90){
  
  res_lls$NRMSE[i] <- nrmse(ion_lls[[i]],ion_miss[[i]],ion[rownames(ion_miss[[i]]),])
  res_lls$MAE.B.A[i] <- mean(abs(sum_lls[[i]]$FC.B.A - sum_lls[[i]]$ori.FC.B.A))
  res_lls$MAE.C.A[i] <- mean(abs(sum_lls[[i]]$FC.C.A - sum_lls[[i]]$ori.FC.C.A))
  res_lls$MAE.D.A[i] <- mean(abs(sum_lls[[i]]$FC.D.A - sum_lls[[i]]$ori.FC.D.A))
  res_lls$NRMSE.B.A[i] <- cal.nrmse(sum_lls[[i]]$ori.FC.B.A,sum_lls[[i]]$FC.B.A)
  res_lls$NRMSE.C.A[i] <- cal.nrmse(sum_lls[[i]]$ori.FC.C.A,sum_lls[[i]]$FC.C.A)
  res_lls$NRMSE.D.A[i] <- cal.nrmse(sum_lls[[i]]$ori.FC.D.A,sum_lls[[i]]$FC.D.A)
  res_lls$TP.B.A[i] <- sum(sum_lls[[i]]$adj.p.B.A < 0.05 & sum_lls[[i]]$species == "ECOLI")
  res_lls$TP.C.A[i] <- sum(sum_lls[[i]]$adj.p.C.A < 0.05 & sum_lls[[i]]$species == "ECOLI")
  res_lls$TP.D.A[i] <- sum(sum_lls[[i]]$adj.p.D.A < 0.05 & sum_lls[[i]]$species == "ECOLI")
  res_lls$FP.B.A[i] <- sum(sum_lls[[i]]$adj.p.B.A < 0.05 & sum_lls[[i]]$species == "HUMAN")
  res_lls$FP.C.A[i] <- sum(sum_lls[[i]]$adj.p.C.A < 0.05 & sum_lls[[i]]$species == "HUMAN")
  res_lls$FP.D.A[i] <- sum(sum_lls[[i]]$adj.p.D.A < 0.05 & sum_lls[[i]]$species == "HUMAN")
  
}

res_lls$FADR.B.A <- res_lls$FP.B.A/(res_lls$TP.B.A+res_lls$FP.B.A)
res_lls$FADR.C.A <- res_lls$FP.C.A/(res_lls$TP.C.A+res_lls$FP.C.A)
res_lls$FADR.D.A <- res_lls$FP.D.A/(res_lls$TP.D.A+res_lls$FP.D.A)

## rf

res_rf <- data.frame(matrix(NA,nrow = 90, ncol = 17))
colnames(res_rf) <- c("time","NRMSE","MAE.B.A","MAE.C.A","MAE.D.A","NRMSE.B.A","NRMSE.C.A","NRMSE.D.A","TP.B.A","TP.C.A","TP.D.A","FP.B.A","FP.C.A","FP.D.A","FADR.B.A","FADR.C.A","FADR.D.A")
rownames(res_rf) <- rownames(time.table)
res_rf$time <- time.table$RF

for(i in 1:90){
  
  res_rf$NRMSE[i] <- nrmse(ion_rf[[i]],ion_miss[[i]],ion[rownames(ion_miss[[i]]),])
  res_rf$MAE.B.A[i] <- mean(abs(sum_rf[[i]]$FC.B.A - sum_rf[[i]]$ori.FC.B.A))
  res_rf$MAE.C.A[i] <- mean(abs(sum_rf[[i]]$FC.C.A - sum_rf[[i]]$ori.FC.C.A))
  res_rf$MAE.D.A[i] <- mean(abs(sum_rf[[i]]$FC.D.A - sum_rf[[i]]$ori.FC.D.A))
  res_rf$NRMSE.B.A[i] <- cal.nrmse(sum_rf[[i]]$ori.FC.B.A,sum_rf[[i]]$FC.B.A)
  res_rf$NRMSE.C.A[i] <- cal.nrmse(sum_rf[[i]]$ori.FC.C.A,sum_rf[[i]]$FC.C.A)
  res_rf$NRMSE.D.A[i] <- cal.nrmse(sum_rf[[i]]$ori.FC.D.A,sum_rf[[i]]$FC.D.A)
  res_rf$TP.B.A[i] <- sum(sum_rf[[i]]$adj.p.B.A < 0.05 & sum_rf[[i]]$species == "ECOLI")
  res_rf$TP.C.A[i] <- sum(sum_rf[[i]]$adj.p.C.A < 0.05 & sum_rf[[i]]$species == "ECOLI")
  res_rf$TP.D.A[i] <- sum(sum_rf[[i]]$adj.p.D.A < 0.05 & sum_rf[[i]]$species == "ECOLI")
  res_rf$FP.B.A[i] <- sum(sum_rf[[i]]$adj.p.B.A < 0.05 & sum_rf[[i]]$species == "HUMAN")
  res_rf$FP.C.A[i] <- sum(sum_rf[[i]]$adj.p.C.A < 0.05 & sum_rf[[i]]$species == "HUMAN")
  res_rf$FP.D.A[i] <- sum(sum_rf[[i]]$adj.p.D.A < 0.05 & sum_rf[[i]]$species == "HUMAN")
  
}

res_rf$FADR.B.A <- res_rf$FP.B.A/(res_rf$TP.B.A+res_rf$FP.B.A)
res_rf$FADR.C.A <- res_rf$FP.C.A/(res_rf$TP.C.A+res_rf$FP.C.A)
res_rf$FADR.D.A <- res_rf$FP.D.A/(res_rf$TP.D.A+res_rf$FP.D.A)

## svd

res_svd <- data.frame(matrix(NA,nrow = 90, ncol = 17))
colnames(res_svd) <- c("time","NRMSE","MAE.B.A","MAE.C.A","MAE.D.A","NRMSE.B.A","NRMSE.C.A","NRMSE.D.A","TP.B.A","TP.C.A","TP.D.A","FP.B.A","FP.C.A","FP.D.A","FADR.B.A","FADR.C.A","FADR.D.A")
rownames(res_svd) <- rownames(time.table)
res_svd$time <- time.table$SVD

for(i in 1:90){
  
  res_svd$NRMSE[i] <- nrmse(ion_svd[[i]],ion_miss[[i]],ion[rownames(ion_miss[[i]]),])
  res_svd$MAE.B.A[i] <- mean(abs(sum_svd[[i]]$FC.B.A - sum_svd[[i]]$ori.FC.B.A))
  res_svd$MAE.C.A[i] <- mean(abs(sum_svd[[i]]$FC.C.A - sum_svd[[i]]$ori.FC.C.A))
  res_svd$MAE.D.A[i] <- mean(abs(sum_svd[[i]]$FC.D.A - sum_svd[[i]]$ori.FC.D.A))
  res_svd$NRMSE.B.A[i] <- cal.nrmse(sum_svd[[i]]$ori.FC.B.A,sum_svd[[i]]$FC.B.A)
  res_svd$NRMSE.C.A[i] <- cal.nrmse(sum_svd[[i]]$ori.FC.C.A,sum_svd[[i]]$FC.C.A)
  res_svd$NRMSE.D.A[i] <- cal.nrmse(sum_svd[[i]]$ori.FC.D.A,sum_svd[[i]]$FC.D.A)
  res_svd$TP.B.A[i] <- sum(sum_svd[[i]]$adj.p.B.A < 0.05 & sum_svd[[i]]$species == "ECOLI")
  res_svd$TP.C.A[i] <- sum(sum_svd[[i]]$adj.p.C.A < 0.05 & sum_svd[[i]]$species == "ECOLI")
  res_svd$TP.D.A[i] <- sum(sum_svd[[i]]$adj.p.D.A < 0.05 & sum_svd[[i]]$species == "ECOLI")
  res_svd$FP.B.A[i] <- sum(sum_svd[[i]]$adj.p.B.A < 0.05 & sum_svd[[i]]$species == "HUMAN")
  res_svd$FP.C.A[i] <- sum(sum_svd[[i]]$adj.p.C.A < 0.05 & sum_svd[[i]]$species == "HUMAN")
  res_svd$FP.D.A[i] <- sum(sum_svd[[i]]$adj.p.D.A < 0.05 & sum_svd[[i]]$species == "HUMAN")
  
}

res_svd$FADR.B.A <- res_svd$FP.B.A/(res_svd$TP.B.A+res_svd$FP.B.A)
res_svd$FADR.C.A <- res_svd$FP.C.A/(res_svd$TP.C.A+res_svd$FP.C.A)
res_svd$FADR.D.A <- res_svd$FP.D.A/(res_svd$TP.D.A+res_svd$FP.D.A)

## bpca

res_bpca <- data.frame(matrix(NA,nrow = 90, ncol = 17))
colnames(res_bpca) <- c("time","NRMSE","MAE.B.A","MAE.C.A","MAE.D.A","NRMSE.B.A","NRMSE.C.A","NRMSE.D.A","TP.B.A","TP.C.A","TP.D.A","FP.B.A","FP.C.A","FP.D.A","FADR.B.A","FADR.C.A","FADR.D.A")
rownames(res_bpca) <- rownames(time.table)
res_bpca$time <- time.table$BPCA

for(i in 1:90){
  
  res_bpca$NRMSE[i] <- nrmse(ion_bpca[[i]],ion_miss[[i]],ion[rownames(ion_miss[[i]]),])
  res_bpca$MAE.B.A[i] <- mean(abs(sum_bpca[[i]]$FC.B.A - sum_bpca[[i]]$ori.FC.B.A))
  res_bpca$MAE.C.A[i] <- mean(abs(sum_bpca[[i]]$FC.C.A - sum_bpca[[i]]$ori.FC.C.A))
  res_bpca$MAE.D.A[i] <- mean(abs(sum_bpca[[i]]$FC.D.A - sum_bpca[[i]]$ori.FC.D.A))
  res_bpca$NRMSE.B.A[i] <- cal.nrmse(sum_bpca[[i]]$ori.FC.B.A,sum_bpca[[i]]$FC.B.A)
  res_bpca$NRMSE.C.A[i] <- cal.nrmse(sum_bpca[[i]]$ori.FC.C.A,sum_bpca[[i]]$FC.C.A)
  res_bpca$NRMSE.D.A[i] <- cal.nrmse(sum_bpca[[i]]$ori.FC.D.A,sum_bpca[[i]]$FC.D.A)
  res_bpca$TP.B.A[i] <- sum(sum_bpca[[i]]$adj.p.B.A < 0.05 & sum_bpca[[i]]$species == "ECOLI")
  res_bpca$TP.C.A[i] <- sum(sum_bpca[[i]]$adj.p.C.A < 0.05 & sum_bpca[[i]]$species == "ECOLI")
  res_bpca$TP.D.A[i] <- sum(sum_bpca[[i]]$adj.p.D.A < 0.05 & sum_bpca[[i]]$species == "ECOLI")
  res_bpca$FP.B.A[i] <- sum(sum_bpca[[i]]$adj.p.B.A < 0.05 & sum_bpca[[i]]$species == "HUMAN")
  res_bpca$FP.C.A[i] <- sum(sum_bpca[[i]]$adj.p.C.A < 0.05 & sum_bpca[[i]]$species == "HUMAN")
  res_bpca$FP.D.A[i] <- sum(sum_bpca[[i]]$adj.p.D.A < 0.05 & sum_bpca[[i]]$species == "HUMAN")
  
}

res_bpca$FADR.B.A <- res_bpca$FP.B.A/(res_bpca$TP.B.A+res_bpca$FP.B.A)
res_bpca$FADR.C.A <- res_bpca$FP.C.A/(res_bpca$TP.C.A+res_bpca$FP.C.A)
res_bpca$FADR.D.A <- res_bpca$FP.D.A/(res_bpca$TP.D.A+res_bpca$FP.D.A)

## export and save results

sum_result <- list("LOD" = res_lod,
                   "ND" = res_nd,
                   "kNN" = res_knn,
                   "LLS" = res_lls,
                   "RF" = res_rf,
                   "SVD" = res_svd,
                   "BPCA" = res_bpca)
library(openxlsx)
write.xlsx(sum_result,file = "20200329_3 proteome summary results.xlsx")

#save(sum_ion,sum_lod,sum_nd,sum_knn,sum_lls,sum_rf,sum_svd,sum_bpca,sum_result,file = "201911_summary and reults.RData")
#load("201911_summary and reults.RData")

#save(sum_ion,sum_lod,sum_nd,sum_knn,sum_lls,sum_rf,sum_svd,sum_bpca,sum_result,file = "202003_summary and reults.RData")
load("202003_summary and reults.RData")

# boxplot -----------------------------------------------------------------

## boxplot of running time

plot_color <- c("red","purple","blue","hotpink","green","deepskyblue","gold")

df_time <- data.frame("Time" = c(sum_result$LOD[,1],sum_result$ND[,1],sum_result$kNN[,1],sum_result$LLS[,1],sum_result$RF[,1],sum_result$SVD[,1],sum_result$BPCA[,1]),
                      "Method" = factor(c(rep("LOD",90),rep("ND",90),rep("kNN",90),rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("LOD","ND","kNN","LLS","RF","SVD","BPCA")),
                      "MV.ratio" = factor(c(rep("10%MV",30),rep("20%MV",30),rep("30%MV",30)),levels = c("10%MV","20%MV","30%MV")),
                      "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",10),rep("50%MNAR",10),rep("80%MNAR",10)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR"))
)

p <- ggplot(df_time,aes(x=Method,y=Time,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("Run time (s)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)


## boxplot of NRMSE

df_nrmse <- data.frame("NRMSE" = c(sum_result$LOD[,2],sum_result$ND[,2],sum_result$kNN[,2],sum_result$LLS[,2],sum_result$RF[,2],sum_result$SVD[,2],sum_result$BPCA[,2]),
                       "Method" = factor(c(rep("LOD",90),rep("ND",90),rep("kNN",90),rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("LOD","ND","kNN","LLS","RF","SVD","BPCA")),
                       "MV.ratio" = factor(c(rep("10%MV",30),rep("20%MV",30),rep("30%MV",30)),levels = c("10%MV","20%MV","30%MV")),
                       "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",10),rep("50%MNAR",10),rep("80%MNAR",10)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR"))
)

p <- ggplot(df_nrmse,aes(x=Method,y=NRMSE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("NRMSE")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)

df_nrmse2 <- df_nrmse[271:630,]
df_nrmse2$Method <- factor(c(rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("kNN","LLS","RF","SVD","BPCA"))

p <- ggplot(df_nrmse2,aes(x=Method,y=NRMSE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color[4:7])
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("NRMSE")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
#p <- p + ylim(0,0.45)
plot(p)



## boxplot of MAE B/A

df_maeB <- data.frame("MAE" = c(sum_result$LOD[,3],sum_result$ND[,3],sum_result$kNN[,3],sum_result$LLS[,3],sum_result$RF[,3],sum_result$SVD[,3],sum_result$BPCA[,3]),
                      "Method" = factor(c(rep("LOD",90),rep("ND",90),rep("kNN",90),rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("LOD","ND","kNN","LLS","RF","SVD","BPCA")),
                      "MV.ratio" = factor(c(rep("10%MV",30),rep("20%MV",30),rep("30%MV",30)),levels = c("10%MV","20%MV","30%MV")),
                      "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",10),rep("50%MNAR",10),rep("80%MNAR",10)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR"))
)

p <- ggplot(df_maeB,aes(x=Method,y=MAE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("MAE (B/A)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)

df_maeB2 <- df_maeB[271:630,]
df_maeB2$Method <- factor(c(rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("kNN","LLS","RF","SVD","BPCA"))

p <- ggplot(df_maeB2,aes(x=Method,y=MAE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color[4:7])
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("MAE (B/A)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
#p <- p + ylim(0,0.18)
plot(p)

## boxplot of MAE C/A

df_maeC <- data.frame("MAE" = c(sum_result$LOD[,4],sum_result$ND[,4],sum_result$kNN[,4],sum_result$LLS[,4],sum_result$RF[,4],sum_result$SVD[,4],sum_result$BPCA[,4]),
                      "Method" = factor(c(rep("LOD",90),rep("ND",90),rep("kNN",90),rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("LOD","ND","kNN","LLS","RF","SVD","BPCA")),
                      "MV.ratio" = factor(c(rep("10%MV",30),rep("20%MV",30),rep("30%MV",30)),levels = c("10%MV","20%MV","30%MV")),
                      "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",10),rep("50%MNAR",10),rep("80%MNAR",10)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR"))
)

p <- ggplot(df_maeC,aes(x=Method,y=MAE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("MAE (C/A)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)

df_maeC2 <- df_maeC[271:630,]
df_maeC2$Method <- factor(c(rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("kNN","LLS","RF","SVD","BPCA"))

p <- ggplot(df_maeC2,aes(x=Method,y=MAE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color[4:7])
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("MAE (C/A)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
#p <- p + ylim(0,0.23)
plot(p)

## boxplot of MAE D/A

df_maeD <- data.frame("MAE" = c(sum_result$LOD[,5],sum_result$ND[,5],sum_result$kNN[,5],sum_result$LLS[,5],sum_result$RF[,5],sum_result$SVD[,5],sum_result$BPCA[,5]),
                      "Method" = factor(c(rep("LOD",90),rep("ND",90),rep("kNN",90),rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("LOD","ND","kNN","LLS","RF","SVD","BPCA")),
                      "MV.ratio" = factor(c(rep("10%MV",30),rep("20%MV",30),rep("30%MV",30)),levels = c("10%MV","20%MV","30%MV")),
                      "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",10),rep("50%MNAR",10),rep("80%MNAR",10)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR"))
)

p <- ggplot(df_maeD,aes(x=Method,y=MAE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("MAE (D/A)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)

df_maeD2 <- df_maeD[271:630,]
df_maeD2$Method <- factor(c(rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("kNN","LLS","RF","SVD","BPCA"))

p <- ggplot(df_maeD2,aes(x=Method,y=MAE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color[4:7])
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("MAE (D/A)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
#p <- p + ylim(0,0.5)
plot(p)


# Boxplot of fold changes -------------------------------------------------

plot_color <- c("red","purple","blue","hotpink","green","deepskyblue","gold")

i = 1

tmp <- do.call(rbind.data.frame, sum_lod[i])
df_lod <- data.frame("value" = c(tmp$FC.B.A,tmp$FC.C.A,tmp$FC.D.A),
                     "Ratio" = c(rep("B/A",nrow(tmp)),rep("C/A",nrow(tmp)),rep("D/A",nrow(tmp))),
                     "Species" = rep(tmp$species,3),
                     "Method" = "LOD")

tmp <- do.call(rbind.data.frame, sum_nd[i])
df_nd <- data.frame("value" = c(tmp$FC.B.A,tmp$FC.C.A,tmp$FC.D.A),
                     "Ratio" = c(rep("B/A",nrow(tmp)),rep("C/A",nrow(tmp)),rep("D/A",nrow(tmp))),
                     "Species" = rep(tmp$species,3),
                     "Method" = "ND")

tmp <- do.call(rbind.data.frame, sum_knn[i])
df_knn <- data.frame("value" = c(tmp$FC.B.A,tmp$FC.C.A,tmp$FC.D.A),
                    "Ratio" = c(rep("B/A",nrow(tmp)),rep("C/A",nrow(tmp)),rep("D/A",nrow(tmp))),
                    "Species" = rep(tmp$species,3),
                    "Method" = "KNN")

tmp <- do.call(rbind.data.frame, sum_lls[i])
df_lls <- data.frame("value" = c(tmp$FC.B.A,tmp$FC.C.A,tmp$FC.D.A),
                     "Ratio" = c(rep("B/A",nrow(tmp)),rep("C/A",nrow(tmp)),rep("D/A",nrow(tmp))),
                     "Species" = rep(tmp$species,3),
                     "Method" = "LLS")

tmp <- do.call(rbind.data.frame, sum_rf[i])
df_rf <- data.frame("value" = c(tmp$FC.B.A,tmp$FC.C.A,tmp$FC.D.A),
                     "Ratio" = c(rep("B/A",nrow(tmp)),rep("C/A",nrow(tmp)),rep("D/A",nrow(tmp))),
                     "Species" = rep(tmp$species,3),
                     "Method" = "RF")

tmp <- do.call(rbind.data.frame, sum_svd[i])
df_svd <- data.frame("value" = c(tmp$FC.B.A,tmp$FC.C.A,tmp$FC.D.A),
                    "Ratio" = c(rep("B/A",nrow(tmp)),rep("C/A",nrow(tmp)),rep("D/A",nrow(tmp))),
                    "Species" = rep(tmp$species,3),
                    "Method" = "SVD")

tmp <- do.call(rbind.data.frame, sum_bpca[i])
df_bpca <- data.frame("value" = c(tmp$FC.B.A,tmp$FC.C.A,tmp$FC.D.A),
                    "Ratio" = c(rep("B/A",nrow(tmp)),rep("C/A",nrow(tmp)),rep("D/A",nrow(tmp))),
                    "Species" = rep(tmp$species,3),
                    "Method" = "BPCA")


df_fc <- rbind(df_lod,df_nd,df_knn,df_lls,df_rf,df_svd,df_bpca)

p <- ggplot(df_fc,aes(x=Method,y=value,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1)# + geom_point(position = position_jitter(width = 0.1))# 
p <- p +  scale_color_manual(values = plot_color)
p <- p + facet_grid(Species~Ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("Protein Ratio")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
p <- p + ylim(0,3)
plot(p)





# ROC plot ----------------------------------------------------------------

#load("201911_summary and reults.RData")
load("201911_imputation data files.RData")
load("202003_summary and reults.RData")


library(pROC)

get_roc <- function(a){
  
  tmp <- a[,15:19]
  tmp$DE <- grepl("ECOLI|YEAST", tmp$species)
  #tmp$DE <- grepl("ECOLI", tmp$species)
  tmp$BG <- grepl("HUMAN", tmp$species)
  
  for(i in 1:3){
    
    tmp2 <- tmp[,c(i,6,7)]
    tmp2 <- tmp2[order(tmp2[,1],decreasing = F),]
    tmp2$TPR <- cumsum(as.numeric(tmp2$DE)) / sum(tmp$DE)
    tmp2$FPR <- cumsum(as.numeric(tmp2$BG)) / sum(tmp$BG)
    tmp <- cbind(tmp,tmp2[rownames(tmp),4:5])
    
  }
  
  colnames(tmp)[8:13] <- paste0(colnames(tmp)[8:13],".",c("B.A","B.A","C.A","C.A","D.A","D.A"))
  
  return(tmp)
}

get_auc <- function(df,rank = "adj.p.B.A"){
  
  return(as.numeric(auc(roc(df[order(df[,grep(rank,colnames(df))],decreasing = F),"DE"],rev(seq(1,nrow(df))),direction = "<"))))
  
}

roc_ion <- get_roc(sum_ion)
roc_lod <- lapply(sum_lod,get_roc)
roc_nd <- lapply(sum_nd,get_roc)
roc_knn <- lapply(sum_knn,get_roc)
roc_lls <- lapply(sum_lls,get_roc)
roc_rf <- lapply(sum_rf,get_roc)
roc_svd <- lapply(sum_svd,get_roc)
roc_bpca <- lapply(sum_bpca,get_roc)


get_auc(roc_ion,rank = "adj.p.B.A")  #0.9146457
get_auc(roc_ion,rank = "adj.p.C.A")  #0.9634159
get_auc(roc_ion,rank = "adj.p.D.A")  #0.9719284

auc_lod <- data.frame("dataset" = rownames(time.table),
                      "AUC.B.A" = unlist(lapply(roc_lod,get_auc, rank="adj.p.B.A")),
                      "AUC.C.A" = unlist(lapply(roc_lod,get_auc, rank="adj.p.C.A")),
                      "AUC.D.A" = unlist(lapply(roc_lod,get_auc, rank="adj.p.D.A"))
)

auc_nd <- data.frame("dataset" = rownames(time.table),
                     "AUC.B.A" = unlist(lapply(roc_nd,get_auc, rank="adj.p.B.A")),
                     "AUC.C.A" = unlist(lapply(roc_nd,get_auc, rank="adj.p.C.A")),
                     "AUC.D.A" = unlist(lapply(roc_nd,get_auc, rank="adj.p.D.A"))
)

auc_knn <- data.frame("dataset" = rownames(time.table),
                      "AUC.B.A" = unlist(lapply(roc_knn,get_auc, rank="adj.p.B.A")),
                      "AUC.C.A" = unlist(lapply(roc_knn,get_auc, rank="adj.p.C.A")),
                      "AUC.D.A" = unlist(lapply(roc_knn,get_auc, rank="adj.p.D.A"))
)

auc_lls <- data.frame("dataset" = rownames(time.table),
                      "AUC.B.A" = unlist(lapply(roc_lls,get_auc, rank="adj.p.B.A")),
                      "AUC.C.A" = unlist(lapply(roc_lls,get_auc, rank="adj.p.C.A")),
                      "AUC.D.A" = unlist(lapply(roc_lls,get_auc, rank="adj.p.D.A"))
)

auc_rf <- data.frame("dataset" = rownames(time.table),
                     "AUC.B.A" = unlist(lapply(roc_rf,get_auc, rank="adj.p.B.A")),
                     "AUC.C.A" = unlist(lapply(roc_rf,get_auc, rank="adj.p.C.A")),
                     "AUC.D.A" = unlist(lapply(roc_rf,get_auc, rank="adj.p.D.A"))
)

auc_svd <- data.frame("dataset" = rownames(time.table),
                      "AUC.B.A" = unlist(lapply(roc_svd,get_auc, rank="adj.p.B.A")),
                      "AUC.C.A" = unlist(lapply(roc_svd,get_auc, rank="adj.p.C.A")),
                      "AUC.D.A" = unlist(lapply(roc_svd,get_auc, rank="adj.p.D.A"))
)

auc_bpca <- data.frame("dataset" = rownames(time.table),
                       "AUC.B.A" = unlist(lapply(roc_bpca,get_auc, rank="adj.p.B.A")),
                       "AUC.C.A" = unlist(lapply(roc_bpca,get_auc, rank="adj.p.C.A")),
                       "AUC.D.A" = unlist(lapply(roc_bpca,get_auc, rank="adj.p.D.A"))
)

auc_result <- list("LOD" = auc_lod,
                   "ND" = auc_nd,
                   "kNN" = auc_knn,
                   "LLS" = auc_lls,
                   "RF" = auc_rf,
                   "SVD" = auc_svd,
                   "BPCA" = auc_bpca)
library(openxlsx)
#write.xlsx(auc_result,file = "201911_3 proteome auc results.xlsx")
write.xlsx(auc_result,file = "20200329_3 proteome auc results.xlsx")

## boxplot of AUCs

plot_color <- c("red","purple","blue","hotpink","green","deepskyblue","gold")

## B vs A

df_auc_B <- data.frame("AUC" = c(auc_result$LOD[,2],auc_result$ND[,2],auc_result$kNN[,2],auc_result$LLS[,2],auc_result$RF[,2],auc_result$SVD[,2],auc_result$BPCA[,2]),
                       "Method" = factor(c(rep("LOD",90),rep("ND",90),rep("kNN",90),rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("LOD","ND","kNN","LLS","RF","SVD","BPCA")),
                       "MV.ratio" = factor(c(rep("10%MV",30),rep("20%MV",30),rep("30%MV",30)),levels = c("10%MV","20%MV","30%MV")),
                       "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",10),rep("50%MNAR",10),rep("80%MNAR",10)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR"))
)

p <- ggplot(df_auc_B,aes(x=Method,y=AUC,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("AUC")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)

df_auc_B2 <- df_auc_B[271:630,]
df_auc_B2$Method <- factor(c(rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("LLS","RF","SVD","BPCA"))

p <- ggplot(df_auc_B2,aes(x=Method,y=AUC,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color[4:7])
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("AUC")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
#p <- p + ylim(0,0.45)
plot(p)

## C vs A

df_auc_C <- data.frame("AUC" = c(auc_result$LOD[,3],auc_result$ND[,3],auc_result$kNN[,3],auc_result$LLS[,3],auc_result$RF[,3],auc_result$SVD[,3],auc_result$BPCA[,3]),
                       "Method" = factor(c(rep("LOD",90),rep("ND",90),rep("kNN",90),rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("LOD","ND","kNN","LLS","RF","SVD","BPCA")),
                       "MV.ratio" = factor(c(rep("10%MV",30),rep("20%MV",30),rep("30%MV",30)),levels = c("10%MV","20%MV","30%MV")),
                       "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",10),rep("50%MNAR",10),rep("80%MNAR",10)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR"))
)

p <- ggplot(df_auc_C,aes(x=Method,y=AUC,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("AUC")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)

df_auc_C2 <- df_auc_C[271:630,]
df_auc_C2$Method <- factor(c(rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("LLS","RF","SVD","BPCA"))

p <- ggplot(df_auc_C2,aes(x=Method,y=AUC,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color[4:7])
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("AUC")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
#p <- p + ylim(0,0.45)
plot(p)

## D vs A

df_auc_D <- data.frame("AUC" = c(auc_result$LOD[,4],auc_result$ND[,4],auc_result$kNN[,4],auc_result$LLS[,4],auc_result$RF[,4],auc_result$SVD[,4],auc_result$BPCA[,4]),
                       "Method" = factor(c(rep("LOD",90),rep("ND",90),rep("kNN",90),rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("LOD","ND","kNN","LLS","RF","SVD","BPCA")),
                       "MV.ratio" = factor(c(rep("10%MV",30),rep("20%MV",30),rep("30%MV",30)),levels = c("10%MV","20%MV","30%MV")),
                       "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",10),rep("50%MNAR",10),rep("80%MNAR",10)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR"))
)

p <- ggplot(df_auc_D,aes(x=Method,y=AUC,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("AUC")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)

df_auc_D2 <- df_auc_D[271:630,]
df_auc_D2$Method <- factor(c(rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("LLS","RF","SVD","BPCA"))

p <- ggplot(df_auc_D2,aes(x=Method,y=AUC,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color[4:7])
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("AUC")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
#p <- p + ylim(0,0.45)
plot(p)


## plot AUC curve

idx <- 1
sum.pro <- as.numeric(unlist(lapply(roc_lod[seq(idx,idx+89,10)],nrow)))

roc.fun <- function(a,b = "TPR",idx = idx){
  
  tmp <- do.call(rbind.data.frame, a[seq(idx,idx+89,10)])
  
  if(b == "TPR"){
    return(c(tmp$TPR.B.A,tmp$TPR.C.A,tmp$TPR.D.A))
  }
  
  if(b == "FPR"){
    return(c(tmp$FPR.B.A,tmp$FPR.C.A,tmp$FPR.D.A))
  }
  
}

df1 <- data.frame("TPR" = roc.fun(roc_lod,b = "TPR", idx = 1),
                  "FPR" = roc.fun(roc_lod,b = "FPR", idx = 1),
                  "group" = c(rep("B.A",sum(sum.pro)),rep("C.A",sum(sum.pro)),rep("D.A",sum(sum.pro))),
                  "MV.ratio" = factor(c(rep("10%MV",sum(sum.pro[1:3])),rep("20%MV",sum(sum.pro[4:6])),rep("30%MV",sum(sum.pro[7:9]))),levels = c("10%MV","20%MV","30%MV")),
                  "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",6689),rep("50%MNAR",6689),rep("80%MNAR",6689)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR")),
                  "method" = "LOD")


df2 <- data.frame("TPR" = roc.fun(roc_nd,b = "TPR", idx = 1),
                  "FPR" = roc.fun(roc_nd,b = "FPR", idx = 1),
                  "group" = c(rep("B.A",sum(sum.pro)),rep("C.A",sum(sum.pro)),rep("D.A",sum(sum.pro))),
                  "MV.ratio" = factor(c(rep("10%MV",sum(sum.pro[1:3])),rep("20%MV",sum(sum.pro[4:6])),rep("30%MV",sum(sum.pro[7:9]))),levels = c("10%MV","20%MV","30%MV")),
                  "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",6689),rep("50%MNAR",6689),rep("80%MNAR",6689)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR")),
                  "method" = "ND")

df3 <- data.frame("TPR" = roc.fun(roc_knn,b = "TPR", idx = 1),
                  "FPR" = roc.fun(roc_knn,b = "FPR", idx = 1),
                  "group" = c(rep("B.A",sum(sum.pro)),rep("C.A",sum(sum.pro)),rep("D.A",sum(sum.pro))),
                  "MV.ratio" = factor(c(rep("10%MV",sum(sum.pro[1:3])),rep("20%MV",sum(sum.pro[4:6])),rep("30%MV",sum(sum.pro[7:9]))),levels = c("10%MV","20%MV","30%MV")),
                  "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",6689),rep("50%MNAR",6689),rep("80%MNAR",6689)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR")),
                  "method" = "kNN")

df4 <- data.frame("TPR" = roc.fun(roc_lls,b = "TPR", idx = 1),
                  "FPR" = roc.fun(roc_lls,b = "FPR", idx = 1),
                  "group" = c(rep("B.A",sum(sum.pro)),rep("C.A",sum(sum.pro)),rep("D.A",sum(sum.pro))),
                  "MV.ratio" = factor(c(rep("10%MV",sum(sum.pro[1:3])),rep("20%MV",sum(sum.pro[4:6])),rep("30%MV",sum(sum.pro[7:9]))),levels = c("10%MV","20%MV","30%MV")),
                  "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",6689),rep("50%MNAR",6689),rep("80%MNAR",6689)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR")),
                  "method" = "LLS")

df5 <- data.frame("TPR" = roc.fun(roc_rf,b = "TPR", idx = 1),
                  "FPR" = roc.fun(roc_rf,b = "FPR", idx = 1),
                  "group" = c(rep("B.A",sum(sum.pro)),rep("C.A",sum(sum.pro)),rep("D.A",sum(sum.pro))),
                  "MV.ratio" = factor(c(rep("10%MV",sum(sum.pro[1:3])),rep("20%MV",sum(sum.pro[4:6])),rep("30%MV",sum(sum.pro[7:9]))),levels = c("10%MV","20%MV","30%MV")),
                  "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",6689),rep("50%MNAR",6689),rep("80%MNAR",6689)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR")),
                  "method" = "RF")

df6 <- data.frame("TPR" = roc.fun(roc_svd,b = "TPR", idx = 1),
                  "FPR" = roc.fun(roc_svd,b = "FPR", idx = 1),
                  "group" = c(rep("B.A",sum(sum.pro)),rep("C.A",sum(sum.pro)),rep("D.A",sum(sum.pro))),
                  "MV.ratio" = factor(c(rep("10%MV",sum(sum.pro[1:3])),rep("20%MV",sum(sum.pro[4:6])),rep("30%MV",sum(sum.pro[7:9]))),levels = c("10%MV","20%MV","30%MV")),
                  "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",6689),rep("50%MNAR",6689),rep("80%MNAR",6689)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR")),
                  "method" = "SVD")

df7 <- data.frame("TPR" = roc.fun(roc_bpca,b = "TPR", idx = 1),
                  "FPR" = roc.fun(roc_bpca,b = "FPR", idx = 1),
                  "group" = c(rep("B.A",sum(sum.pro)),rep("C.A",sum(sum.pro)),rep("D.A",sum(sum.pro))),
                  "MV.ratio" = factor(c(rep("10%MV",sum(sum.pro[1:3])),rep("20%MV",sum(sum.pro[4:6])),rep("30%MV",sum(sum.pro[7:9]))),levels = c("10%MV","20%MV","30%MV")),
                  "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",6689),rep("50%MNAR",6689),rep("80%MNAR",6689)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR")),
                  "method" = "BPCA")

df.all <- rbind(df1,df2,df3,df4,df5,df6,df7)

library(ggplot2)

p <- ggplot(df.all[df.all$group == "D.A",], aes(FPR, TPR, col = method))
p <- p + geom_line() + scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + theme_bw(base_size = 14) + ggtitle("ROC-curve (D/A)") + xlim(0, 0.1)
plot(p)





# summary table -----------------------------------------------------------

res1 <- as.data.frame(matrix(data = NA, ncol = 8, nrow = 7))
for(i in 1:7){
  res1[i,1] <- mean(sum_result[[i]]$NRMSE)
  res1[i,2] <- sd(sum_result[[i]]$NRMSE)
  res1[i,3] <- mean(sum_result[[i]]$MAE.B.A)
  res1[i,4] <- sd(sum_result[[i]]$MAE.B.A)
  res1[i,5] <- mean(sum_result[[i]]$MAE.C.A)
  res1[i,6] <- sd(sum_result[[i]]$MAE.C.A)
  res1[i,7] <- mean(sum_result[[i]]$MAE.D.A)
  res1[i,8] <- sd(sum_result[[i]]$MAE.D.A)
}

res2<- as.data.frame(matrix(data = NA, ncol = 6, nrow = 7))

for(i in 1:7){
  res2[i,1] <- mean(auc_result[[i]]$AUC.B.A)
  res2[i,2] <- sd(auc_result[[i]]$AUC.B.A)
  res2[i,3] <- mean(auc_result[[i]]$AUC.C.A)
  res2[i,4] <- sd(auc_result[[i]]$AUC.C.A)
  res2[i,5] <- mean(auc_result[[i]]$AUC.D.A)
  res2[i,6] <- sd(auc_result[[i]]$AUC.D.A)
}

res3 <- as.data.frame(matrix(data = NA, ncol = 2, nrow = 7))

for(i in 1:7){
  res3[i,1] <- mean(as.numeric(c(sum_result[[i]]$MAE.B.A,sum_result[[i]]$MAE.C.A,sum_result[[i]]$MAE.D.A)))
  res3[i,2] <- sd(as.numeric(c(sum_result[[i]]$MAE.B.A,sum_result[[i]]$MAE.C.A,sum_result[[i]]$MAE.D.A)))
}

res4 <- as.data.frame(matrix(data = NA, ncol = 2, nrow = 7))

for(i in 1:7){
  res4[i,1] <- mean(as.numeric(c(auc_result[[i]]$AUC.B.A,auc_result[[i]]$AUC.C.A,auc_result[[i]]$AUC.D.A)))
  res4[i,2] <- sd(as.numeric(c(auc_result[[i]]$AUC.B.A,auc_result[[i]]$AUC.C.A,auc_result[[i]]$AUC.D.A)))
}


res1 <- as.data.frame(matrix(data = NA, ncol = 8, nrow = 7))
for(i in 1:7){
  res1[i,1] <- mean(sum_result[[i]]$NRMSE)
  res1[i,2] <- sd(sum_result[[i]]$NRMSE)
  res1[i,3] <- mean(sum_result[[i]]$NRMSE.B.A)
  res1[i,4] <- sd(sum_result[[i]]$NRMSE.B.A)
  res1[i,5] <- mean(sum_result[[i]]$NRMSE.C.A)
  res1[i,6] <- sd(sum_result[[i]]$NRMSE.C.A)
  res1[i,7] <- mean(sum_result[[i]]$NRMSE.D.A)
  res1[i,8] <- sd(sum_result[[i]]$NRMSE.D.A)
}

res3 <- as.data.frame(matrix(data = NA, ncol = 2, nrow = 7))

for(i in 1:7){
  res3[i,1] <- mean(as.numeric(c(sum_result[[i]]$NRMSE.B.A,sum_result[[i]]$NRMSE.C.A,sum_result[[i]]$NRMSE.D.A)))
  res3[i,2] <- sd(as.numeric(c(sum_result[[i]]$NRMSE.B.A,sum_result[[i]]$NRMSE.C.A,sum_result[[i]]$NRMSE.D.A)))
}


res1 <- as.data.frame(matrix(data = NA, ncol = 12, nrow = 7))
for(i in 1:7){
  res1[i,1] <- mean(sum_result[[i]]$TP.B.A)
  res1[i,2] <- sd(sum_result[[i]]$TP.B.A)
  res1[i,3] <- mean(sum_result[[i]]$TP.C.A)
  res1[i,4] <- sd(sum_result[[i]]$TP.C.A)
  res1[i,5] <- mean(sum_result[[i]]$TP.D.A)
  res1[i,6] <- sd(sum_result[[i]]$TP.D.A)
  res1[i,7] <- mean(sum_result[[i]]$FP.B.A)
  res1[i,8] <- sd(sum_result[[i]]$FP.B.A)
  res1[i,9] <- mean(sum_result[[i]]$FP.C.A)
  res1[i,10] <- sd(sum_result[[i]]$FP.C.A)
  res1[i,11] <- mean(sum_result[[i]]$FP.D.A)
  res1[i,12] <- sd(sum_result[[i]]$FP.D.A)
  
}

res3 <- as.data.frame(matrix(data = NA, ncol = 4, nrow = 7))

for(i in 1:7){
  res3[i,1] <- mean(as.numeric(c(sum_result[[i]]$TP.B.A,sum_result[[i]]$TP.C.A,sum_result[[i]]$TP.D.A)))
  res3[i,2] <- sd(as.numeric(c(sum_result[[i]]$TP.B.A,sum_result[[i]]$TP.C.A,sum_result[[i]]$TP.D.A)))
  res3[i,3] <- mean(as.numeric(c(sum_result[[i]]$FP.B.A,sum_result[[i]]$FP.C.A,sum_result[[i]]$FP.D.A)))
  res3[i,4] <- sd(as.numeric(c(sum_result[[i]]$FP.B.A,sum_result[[i]]$FP.C.A,sum_result[[i]]$FP.D.A)))
}


res_lod[res_lod == "NaN"] <- NA
sum_result$LOD <- res_lod

res1 <- as.data.frame(matrix(data = NA, ncol = 6, nrow = 7))
for(i in 1:7){
  res1[i,1] <- mean(sum_result[[i]]$FADR.B.A,na.rm = T)
  res1[i,2] <- sd(sum_result[[i]]$FADR.B.A,na.rm = T)
  res1[i,3] <- mean(sum_result[[i]]$FADR.C.A)
  res1[i,4] <- sd(sum_result[[i]]$FADR.C.A)
  res1[i,5] <- mean(sum_result[[i]]$FADR.D.A)
  res1[i,6] <- sd(sum_result[[i]]$FADR.D.A)
  
}

res4 <- as.data.frame(matrix(data = NA, ncol = 12, nrow = 9))
tmp <- sum_result$BPCA
for(i in 1:9){
  res4[i,1] <- mean(tmp$TP.B.A[(10*(i-1)+1):(10*i)],na.rm = T)
  res4[i,2] <- sd(tmp$TP.B.A[(10*(i-1)+1):(10*i)],na.rm = T)
  res4[i,3] <- mean(tmp$TP.C.A[(10*(i-1)+1):(10*i)],na.rm = T)
  res4[i,4] <- sd(tmp$TP.C.A[(10*(i-1)+1):(10*i)],na.rm = T)
  res4[i,5] <- mean(tmp$TP.D.A[(10*(i-1)+1):(10*i)],na.rm = T)
  res4[i,6] <- sd(tmp$TP.D.A[(10*(i-1)+1):(10*i)],na.rm = T)
  res4[i,7] <- mean(tmp$FADR.B.A[(10*(i-1)+1):(10*i)],na.rm = T)
  res4[i,8] <- sd(tmp$FADR.B.A[(10*(i-1)+1):(10*i)],na.rm = T)
  res4[i,9] <- mean(tmp$FADR.C.A[(10*(i-1)+1):(10*i)],na.rm = T)
  res4[i,10] <- sd(tmp$FADR.C.A[(10*(i-1)+1):(10*i)],na.rm = T)
  res4[i,11] <- mean(tmp$FADR.D.A[(10*(i-1)+1):(10*i)],na.rm = T)
  res4[i,12] <- sd(tmp$FADR.D.A[(10*(i-1)+1):(10*i)],na.rm = T)
}








