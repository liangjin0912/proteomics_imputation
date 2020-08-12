setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# load data and library ---------------------------------------------------

load("201912_imputation data files_cell.RData")
source("functions_cell lineage.R")

library(dplyr)
library(missForest)
library(purrr)
library(ggplot2)

# generate summary tables -------------------------------------------------

sum_ori <- SumTable(ori2)

sum_lod <- lapply(cell_lod,SumTable)
names(sum_lod) <- names(cell_lod)

sum_nd <- lapply(cell_nd,SumTable)
names(sum_nd) <- names(cell_nd)

sum_knn <- lapply(cell_knn,SumTable)
names(sum_knn) <- names(cell_knn)

sum_lls <- lapply(cell_lls,SumTable)
names(sum_lls) <- names(cell_lls)

sum_rf <- lapply(cell_rf,SumTable)
names(sum_rf) <- names(cell_rf)

sum_svd <- lapply(cell_svd,SumTable)
names(sum_svd) <- names(cell_svd)

sum_bpca <- lapply(cell_bpca,SumTable)
names(sum_bpca) <- names(cell_bpca)

# generate result tables --------------------------------------------------

## lod

res_lod <- data.frame(matrix(NA,nrow = 90, ncol = 10))
colnames(res_lod) <- c("time","NRMSE","MAE.B","MAE.T4","MAE.T8","MAE.MO","NRMSE.B","NRMSE.T4","NRMSE.T8","NRMSE.MO")
rownames(res_lod) <- rownames(time.table)
res_lod$time <- time.table$LOD

for(i in 1:90){
  
  res_lod$NRMSE[i] <- nrmse(cell_lod[[i]],cell_miss[[i]],ori2[rownames(cell_miss[[i]]),])
  res_lod$MAE.B[i] <- mean(abs(sum_lod[[i]]$FC.B - sum_lod[[i]]$ori.FC.B))
  res_lod$MAE.T4[i] <- mean(abs(sum_lod[[i]]$FC.T4 - sum_lod[[i]]$ori.FC.T4))
  res_lod$MAE.T8[i] <- mean(abs(sum_lod[[i]]$FC.T8 - sum_lod[[i]]$ori.FC.T8))
  res_lod$MAE.MO[i] <- mean(abs(sum_lod[[i]]$FC.MO - sum_lod[[i]]$ori.FC_MO))
  res_lod$NRMSE.B[i] <- cal.nrmse(sum_lod[[i]]$FC.B,sum_lod[[i]]$ori.FC.B)
  res_lod$NRMSE.T4[i] <- cal.nrmse(sum_lod[[i]]$FC.T4,sum_lod[[i]]$ori.FC.T4)
  res_lod$NRMSE.T8[i] <- cal.nrmse(sum_lod[[i]]$FC.T8,sum_lod[[i]]$ori.FC.T8)
  res_lod$NRMSE.MO[i] <- cal.nrmse(sum_lod[[i]]$FC.MO,sum_lod[[i]]$ori.FC_MO)
  
}

## nd

res_nd <- data.frame(matrix(NA,nrow = 90, ncol = 10))
colnames(res_nd) <- c("time","NRMSE","MAE.B","MAE.T4","MAE.T8","MAE.MO","NRMSE.B","NRMSE.T4","NRMSE.T8","NRMSE.MO")
rownames(res_nd) <- rownames(time.table)
res_nd$time <- time.table$ND

for(i in 1:90){
  
  res_nd$NRMSE[i] <- nrmse(cell_nd[[i]],cell_miss[[i]],ori2[rownames(cell_miss[[i]]),])
  res_nd$MAE.B[i] <- mean(abs(sum_nd[[i]]$FC.B - sum_nd[[i]]$ori.FC.B))
  res_nd$MAE.T4[i] <- mean(abs(sum_nd[[i]]$FC.T4 - sum_nd[[i]]$ori.FC.T4))
  res_nd$MAE.T8[i] <- mean(abs(sum_nd[[i]]$FC.T8 - sum_nd[[i]]$ori.FC.T8))
  res_nd$MAE.MO[i] <- mean(abs(sum_nd[[i]]$FC.MO - sum_nd[[i]]$ori.FC_MO))
  res_nd$NRMSE.B[i] <- cal.nrmse(sum_nd[[i]]$FC.B,sum_nd[[i]]$ori.FC.B)
  res_nd$NRMSE.T4[i] <- cal.nrmse(sum_nd[[i]]$FC.T4,sum_nd[[i]]$ori.FC.T4)
  res_nd$NRMSE.T8[i] <- cal.nrmse(sum_nd[[i]]$FC.T8,sum_nd[[i]]$ori.FC.T8)
  res_nd$NRMSE.MO[i] <- cal.nrmse(sum_nd[[i]]$FC.MO,sum_nd[[i]]$ori.FC_MO)
  
}

## knn

res_knn <- data.frame(matrix(NA,nrow = 90, ncol = 10))
colnames(res_knn) <- c("time","NRMSE","MAE.B","MAE.T4","MAE.T8","MAE.MO","NRMSE.B","NRMSE.T4","NRMSE.T8","NRMSE.MO")
rownames(res_knn) <- rownames(time.table)
res_knn$time <- time.table$kNN

for(i in 1:90){
  
  res_knn$NRMSE[i] <- nrmse(cell_knn[[i]],cell_miss[[i]],ori2[rownames(cell_miss[[i]]),])
  res_knn$MAE.B[i] <- mean(abs(sum_knn[[i]]$FC.B - sum_knn[[i]]$ori.FC.B))
  res_knn$MAE.T4[i] <- mean(abs(sum_knn[[i]]$FC.T4 - sum_knn[[i]]$ori.FC.T4))
  res_knn$MAE.T8[i] <- mean(abs(sum_knn[[i]]$FC.T8 - sum_knn[[i]]$ori.FC.T8))
  res_knn$MAE.MO[i] <- mean(abs(sum_knn[[i]]$FC.MO - sum_knn[[i]]$ori.FC_MO))
  res_knn$NRMSE.B[i] <- cal.nrmse(sum_knn[[i]]$FC.B,sum_knn[[i]]$ori.FC.B)
  res_knn$NRMSE.T4[i] <- cal.nrmse(sum_knn[[i]]$FC.T4,sum_knn[[i]]$ori.FC.T4)
  res_knn$NRMSE.T8[i] <- cal.nrmse(sum_knn[[i]]$FC.T8,sum_knn[[i]]$ori.FC.T8)
  res_knn$NRMSE.MO[i] <- cal.nrmse(sum_knn[[i]]$FC.MO,sum_knn[[i]]$ori.FC_MO)
  
}

## lls

res_lls <- data.frame(matrix(NA,nrow = 90, ncol = 10))
colnames(res_lls) <- c("time","NRMSE","MAE.B","MAE.T4","MAE.T8","MAE.MO","NRMSE.B","NRMSE.T4","NRMSE.T8","NRMSE.MO")
rownames(res_lls) <- rownames(time.table)
res_lls$time <- time.table$LLS

for(i in 1:90){
  
  res_lls$NRMSE[i] <- nrmse(cell_lls[[i]],cell_miss[[i]],ori2[rownames(cell_miss[[i]]),])
  res_lls$MAE.B[i] <- mean(abs(sum_lls[[i]]$FC.B - sum_lls[[i]]$ori.FC.B))
  res_lls$MAE.T4[i] <- mean(abs(sum_lls[[i]]$FC.T4 - sum_lls[[i]]$ori.FC.T4))
  res_lls$MAE.T8[i] <- mean(abs(sum_lls[[i]]$FC.T8 - sum_lls[[i]]$ori.FC.T8))
  res_lls$MAE.MO[i] <- mean(abs(sum_lls[[i]]$FC.MO - sum_lls[[i]]$ori.FC_MO))
  res_lls$NRMSE.B[i] <- cal.nrmse(sum_lls[[i]]$FC.B,sum_lls[[i]]$ori.FC.B)
  res_lls$NRMSE.T4[i] <- cal.nrmse(sum_lls[[i]]$FC.T4,sum_lls[[i]]$ori.FC.T4)
  res_lls$NRMSE.T8[i] <- cal.nrmse(sum_lls[[i]]$FC.T8,sum_lls[[i]]$ori.FC.T8)
  res_lls$NRMSE.MO[i] <- cal.nrmse(sum_lls[[i]]$FC.MO,sum_lls[[i]]$ori.FC_MO)
  
}

## rf

res_rf <- data.frame(matrix(NA,nrow = 90, ncol = 10))
colnames(res_rf) <- c("time","NRMSE","MAE.B","MAE.T4","MAE.T8","MAE.MO","NRMSE.B","NRMSE.T4","NRMSE.T8","NRMSE.MO")
rownames(res_rf) <- rownames(time.table)
res_rf$time <- time.table$RF

for(i in 1:90){
  
  res_rf$NRMSE[i] <- nrmse(cell_rf[[i]],cell_miss[[i]],ori2[rownames(cell_miss[[i]]),])
  res_rf$MAE.B[i] <- mean(abs(sum_rf[[i]]$FC.B - sum_rf[[i]]$ori.FC.B))
  res_rf$MAE.T4[i] <- mean(abs(sum_rf[[i]]$FC.T4 - sum_rf[[i]]$ori.FC.T4))
  res_rf$MAE.T8[i] <- mean(abs(sum_rf[[i]]$FC.T8 - sum_rf[[i]]$ori.FC.T8))
  res_rf$MAE.MO[i] <- mean(abs(sum_rf[[i]]$FC.MO - sum_rf[[i]]$ori.FC_MO))
  res_rf$NRMSE.B[i] <- cal.nrmse(sum_rf[[i]]$FC.B,sum_rf[[i]]$ori.FC.B)
  res_rf$NRMSE.T4[i] <- cal.nrmse(sum_rf[[i]]$FC.T4,sum_rf[[i]]$ori.FC.T4)
  res_rf$NRMSE.T8[i] <- cal.nrmse(sum_rf[[i]]$FC.T8,sum_rf[[i]]$ori.FC.T8)
  res_rf$NRMSE.MO[i] <- cal.nrmse(sum_rf[[i]]$FC.MO,sum_rf[[i]]$ori.FC_MO)
  
}

## svd

res_svd <- data.frame(matrix(NA,nrow = 90, ncol = 10))
colnames(res_svd) <- c("time","NRMSE","MAE.B","MAE.T4","MAE.T8","MAE.MO","NRMSE.B","NRMSE.T4","NRMSE.T8","NRMSE.MO")
rownames(res_svd) <- rownames(time.table)
res_svd$time <- time.table$SVD

for(i in 1:90){
  
  res_svd$NRMSE[i] <- nrmse(cell_svd[[i]],cell_miss[[i]],ori2[rownames(cell_miss[[i]]),])
  res_svd$MAE.B[i] <- mean(abs(sum_svd[[i]]$FC.B - sum_svd[[i]]$ori.FC.B))
  res_svd$MAE.T4[i] <- mean(abs(sum_svd[[i]]$FC.T4 - sum_svd[[i]]$ori.FC.T4))
  res_svd$MAE.T8[i] <- mean(abs(sum_svd[[i]]$FC.T8 - sum_svd[[i]]$ori.FC.T8))
  res_svd$MAE.MO[i] <- mean(abs(sum_svd[[i]]$FC.MO - sum_svd[[i]]$ori.FC_MO))
  res_svd$NRMSE.B[i] <- cal.nrmse(sum_svd[[i]]$FC.B,sum_svd[[i]]$ori.FC.B)
  res_svd$NRMSE.T4[i] <- cal.nrmse(sum_svd[[i]]$FC.T4,sum_svd[[i]]$ori.FC.T4)
  res_svd$NRMSE.T8[i] <- cal.nrmse(sum_svd[[i]]$FC.T8,sum_svd[[i]]$ori.FC.T8)
  res_svd$NRMSE.MO[i] <- cal.nrmse(sum_svd[[i]]$FC.MO,sum_svd[[i]]$ori.FC_MO)
  
}

## bpca

res_bpca <- data.frame(matrix(NA,nrow = 90, ncol = 10))
colnames(res_bpca) <- c("time","NRMSE","MAE.B","MAE.T4","MAE.T8","MAE.MO","NRMSE.B","NRMSE.T4","NRMSE.T8","NRMSE.MO")
rownames(res_bpca) <- rownames(time.table)
res_bpca$time <- time.table$BPCA

for(i in 1:90){
  
  res_bpca$NRMSE[i] <- nrmse(cell_bpca[[i]],cell_miss[[i]],ori2[rownames(cell_miss[[i]]),])
  res_bpca$MAE.B[i] <- mean(abs(sum_bpca[[i]]$FC.B - sum_bpca[[i]]$ori.FC.B))
  res_bpca$MAE.T4[i] <- mean(abs(sum_bpca[[i]]$FC.T4 - sum_bpca[[i]]$ori.FC.T4))
  res_bpca$MAE.T8[i] <- mean(abs(sum_bpca[[i]]$FC.T8 - sum_bpca[[i]]$ori.FC.T8))
  res_bpca$MAE.MO[i] <- mean(abs(sum_bpca[[i]]$FC.MO - sum_bpca[[i]]$ori.FC_MO))
  res_bpca$NRMSE.B[i] <- cal.nrmse(sum_bpca[[i]]$FC.B,sum_bpca[[i]]$ori.FC.B)
  res_bpca$NRMSE.T4[i] <- cal.nrmse(sum_bpca[[i]]$FC.T4,sum_bpca[[i]]$ori.FC.T4)
  res_bpca$NRMSE.T8[i] <- cal.nrmse(sum_bpca[[i]]$FC.T8,sum_bpca[[i]]$ori.FC.T8)
  res_bpca$NRMSE.MO[i] <- cal.nrmse(sum_bpca[[i]]$FC.MO,sum_bpca[[i]]$ori.FC_MO)
  
}

## export and save results

sum_result <- list("LOD" = res_lod,
                   "ND" = res_nd,
                   "kNN" = res_knn,
                   "LLS" = res_lls,
                   "RF" = res_rf,
                   "SVD" = res_svd,
                   "BPCA" = res_bpca)
library(openxlsx)
write.xlsx(sum_result,file = "202003_cell lineage summary results.xlsx")

#save(sum_ori,sum_lod,sum_nd,sum_knn,sum_lls,sum_rf,sum_svd,sum_bpca,sum_result,file = "202001_summary and reults.RData")
#load("202001_summary and reults.RData")
#save(sum_ori,sum_lod,sum_nd,sum_knn,sum_lls,sum_rf,sum_svd,sum_bpca,sum_result,file = "202003_summary and reults.RData")
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
p <- p + xlab("Imputation methods") + ylab("MAE (B cells)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)

df_maeB2 <- df_maeB[271:630,]
df_maeB2$Method <- factor(c(rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("kNN","LLS","RF","SVD","BPCA"))

p <- ggplot(df_maeB2,aes(x=Method,y=MAE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color[4:7])
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("MAE (B cells)")
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
p <- p + xlab("Imputation methods") + ylab("MAE (T4 cells)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)

df_maeC2 <- df_maeC[271:630,]
df_maeC2$Method <- factor(c(rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("kNN","LLS","RF","SVD","BPCA"))

p <- ggplot(df_maeC2,aes(x=Method,y=MAE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color[4:7])
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("MAE (T4 cells)")
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
p <- p + xlab("Imputation methods") + ylab("MAE (T8 cells)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)

df_maeD2 <- df_maeD[271:630,]
df_maeD2$Method <- factor(c(rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("kNN","LLS","RF","SVD","BPCA"))

p <- ggplot(df_maeD2,aes(x=Method,y=MAE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color[4:7])
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("MAE (T8 cells)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
#p <- p + ylim(0,0.5)
plot(p)


## boxplot of MAE MO

df_maeE <- data.frame("MAE" = c(sum_result$LOD[,6],sum_result$ND[,6],sum_result$kNN[,6],sum_result$LLS[,6],sum_result$RF[,6],sum_result$SVD[,6],sum_result$BPCA[,6]),
                      "Method" = factor(c(rep("LOD",90),rep("ND",90),rep("kNN",90),rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("LOD","ND","kNN","LLS","RF","SVD","BPCA")),
                      "MV.ratio" = factor(c(rep("10%MV",30),rep("20%MV",30),rep("30%MV",30)),levels = c("10%MV","20%MV","30%MV")),
                      "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",10),rep("50%MNAR",10),rep("80%MNAR",10)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR"))
)

p <- ggplot(df_maeE,aes(x=Method,y=MAE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("MAE (Monocytes)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)

df_maeE2 <- df_maeE[271:630,]
df_maeE2$Method <- factor(c(rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("kNN","LLS","RF","SVD","BPCA"))

p <- ggplot(df_maeE2,aes(x=Method,y=MAE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color[4:7])
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("MAE (Monocytes)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
#p <- p + ylim(0,0.5)
plot(p)


# ROC plot ----------------------------------------------------------------

load("202001_summary and reults.RData")
load("201912_imputation data files_cell.RData")


library(pROC)

get_roc <- function(a){
  
  tmp <- a[,15:19]
  tmp$DE <- grepl("ECOLI|YEAST", tmp$species)
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
write.xlsx(auc_result,file = "201911_3 proteome auc results.xlsx")

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



# summary table -----------------------------------------------------------

res1 <- as.data.frame(matrix(data = NA, ncol = 10, nrow = 7))
for(i in 1:7){
  res1[i,1] <- mean(sum_result[[i]]$NRMSE)
  res1[i,2] <- sd(sum_result[[i]]$NRMSE)
  res1[i,3] <- mean(sum_result[[i]]$MAE.B)
  res1[i,4] <- sd(sum_result[[i]]$MAE.B)
  res1[i,5] <- mean(sum_result[[i]]$MAE.T4)
  res1[i,6] <- sd(sum_result[[i]]$MAE.T4)
  res1[i,7] <- mean(sum_result[[i]]$MAE.T8)
  res1[i,8] <- sd(sum_result[[i]]$MAE.T8)
  res1[i,9] <- mean(sum_result[[i]]$MAE.MO)
  res1[i,10] <- sd(sum_result[[i]]$MAE.MO)
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
  res3[i,1] <- mean(as.numeric(c(sum_result[[i]]$MAE.B,sum_result[[i]]$MAE.T4,sum_result[[i]]$MAE.T8,sum_result[[i]]$MAE.MO)))
  res3[i,2] <- sd(as.numeric(c(sum_result[[i]]$MAE.B,sum_result[[i]]$MAE.T4,sum_result[[i]]$MAE.T8,sum_result[[i]]$MAE.MO)))
}


res1 <- as.data.frame(matrix(data = NA, ncol = 10, nrow = 7))
for(i in 1:7){
  res1[i,1] <- mean(sum_result[[i]]$NRMSE)
  res1[i,2] <- sd(sum_result[[i]]$NRMSE)
  res1[i,3] <- mean(sum_result[[i]]$NRMSE.B)
  res1[i,4] <- sd(sum_result[[i]]$NRMSE.B)
  res1[i,5] <- mean(sum_result[[i]]$NRMSE.T4)
  res1[i,6] <- sd(sum_result[[i]]$NRMSE.T4)
  res1[i,7] <- mean(sum_result[[i]]$NRMSE.T8)
  res1[i,8] <- sd(sum_result[[i]]$NRMSE.T8)
  res1[i,9] <- mean(sum_result[[i]]$NRMSE.MO)
  res1[i,10] <- sd(sum_result[[i]]$NRMSE.MO)
}

res3 <- as.data.frame(matrix(data = NA, ncol = 2, nrow = 7))

for(i in 1:7){
  res3[i,1] <- mean(as.numeric(c(sum_result[[i]]$NRMSE.B,sum_result[[i]]$NRMSE.T4,sum_result[[i]]$NRMSE.T8,sum_result[[i]]$NRMSE.MO)))
  res3[i,2] <- sd(as.numeric(c(sum_result[[i]]$NRMSE.B,sum_result[[i]]$NRMSE.T4,sum_result[[i]]$NRMSE.T8,sum_result[[i]]$NRMSE.MO)))
}