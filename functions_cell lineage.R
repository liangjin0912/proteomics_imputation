# functions ---------------------------------------------------------------

## calculate NRMSE

cal.nrmse <- function(sim,obs){
  
  sqrt(sum((sim-obs)^2)/length(sim))/mean(c(sim,obs))
  
}

## input dataframe has to be formatted as ionstar data

SumTable <- function(res){
  
  res <-res[!grepl("TRYP_PIG",rownames(res)),]
  
  a <- 2^res[,1:4]
  b <- 2^res[,5:8]
  c <- 2^res[,9:12]
  d <- 2^res[,13:16]
  e <- 2^res[,17:20]
  f <- 2^res[,21:24]
  g <- 2^res[,25:28]
  h <- 2^res[,29:31]
  
  tmp1 <- data.frame("Ave.B.steady" = apply(a,1,mean),
                     "Ave.T4.steady" = apply(b,1,mean),
                     "Ave.T8.steady" = apply(c,1,mean),
                     "Ave.MO.steady" = apply(d,1,mean),
                     "Ave.B.activated" = apply(e,1,mean),
                     "Ave.T4.activated" = apply(f,1,mean),
                     "Ave.T8.activated" = apply(g,1,mean),
                     "Ave.MO.activated" = apply(h,1,mean)
                     )
  
  tmp2 <- data.frame("Stdev.B.steady" = apply(a,1,sd),
                     "Stdev.T4.steady" = apply(b,1,sd),
                     "Stdev.T8.steady" = apply(c,1,sd),
                     "Stdev.MO.steady" = apply(d,1,sd),
                     "Stdev.B.activated" = apply(e,1,sd),
                     "Stdev.T4.activated" = apply(f,1,sd),
                     "Stdev.T8.activated" = apply(g,1,sd),
                     "Stdev.MO.activated" = apply(h,1,sd)
                     )
  
  tmp3 <- tmp2/tmp1
  colnames(tmp3) <- c("CV.B.steady","CV.T4.steady","CV.T8.steady","CV.MO.steady","CV.B.activated","CV.T4.activated","CV.T8.activated","CV.MO.activated")
  
  tmp4 <- data.frame("FC.B" = tmp1$Ave.B.activated / tmp1$Ave.B.steady,
                     "FC.T4" = tmp1$Ave.T4.activated / tmp1$Ave.T4.steady,
                     "FC.T8" = tmp1$Ave.T8.activated / tmp1$Ave.T8.steady,
                     "FC.MO" = tmp1$Ave.MO.activated / tmp1$Ave.MO.steady)
  
  a <- res[,1:4]
  b <- res[,5:8]
  c <- res[,9:12]
  d <- res[,13:16]
  e <- res[,17:20]
  f <- res[,21:24]
  g <- res[,25:28]
  h <- res[,29:31]
  
  tmp5 <- data.frame("p.B" = apply(cbind(a,e),1,getp),
                     "p.T4" = apply(cbind(b,f),1,getp),
                     "p.T8" = apply(cbind(c,g),1,getp),
                     "p.MO" = apply(cbind(d,h),1,getp)
                     )

  rownames(tmp4) <- rownames(res)
  
  tmp6 <- data.frame("adj.p.B" = p.adjust(tmp5$p.B,method = "BH"),
                     "adj.p.T4" = p.adjust(tmp5$p.T4,method = "BH"),
                     "adj.p.T8" = p.adjust(tmp5$p.T8,method = "BH"),
                     "adj.p.MO" = p.adjust(tmp5$p.MO,method = "BH"),stringsAsFactors = F)
  
  tmp7 <- data.frame("ori.FC.B" = apply(2^(ori2[,17:20]),1,mean)/apply(2^(ori2[,1:4]),1,mean),
                     "ori.FC.T4" = apply(2^(ori2[,21:24]),1,mean)/apply(2^(ori2[,5:8]),1,mean),
                     "ori.FC.T8" = apply(2^(ori2[,25:28]),1,mean)/apply(2^(ori2[,9:12]),1,mean),
                     "ori.FC_MO" = apply(2^(ori2[,29:31]),1,mean)/apply(2^(ori2[,13:16]),1,mean)
                     )
  rownames(tmp7) <- rownames(ori2)
  
  return(cbind(tmp1,tmp3,tmp4,tmp5,tmp6,tmp7[rownames(res),]))
  
}

## get original fold change from a vector

orifc <- function(x, human = 1, ecoli = 1, yeast = 1){
  
  x[x == "HUMAN"] <- human
  x[x == "ECOLI"] <- ecoli
  x[x == "YEAST"] <- yeast
  
  return(as.numeric(x))
}

## get p-value

getp <- function(res2){
  if(length(res2) < 8){
    myTtest <- try(t.test(res2[1:4],res2[5:7],var.equal = T,conf.level=0.95))
  }else{
    myTtest <- try(t.test(res2[1:4],res2[5:8],var.equal = T,conf.level=0.95))
  }
  
  if (inherits(myTtest, "try-error"))
  {
    return(1)
  }else{
    return(myTtest$p.value)
  }
}

## summarize TP, FP and FADR from sumtable function output
## bug fixed 12/04/18

fadr <- function(res3){
  
  tmp1 <- matrix(nrow = 4,ncol = 4)
  tmp1[1,1] <- sum(grepl("ECOLI",rownames(res3[res3$FC.B.A >= 1.4 & res3$p.B.A < 0.05,])))
  tmp1[1,2] <- sum(grepl("HUMAN",rownames(res3[(res3$FC.B.A >= 1.4 & res3$p.B.A < 0.05)|(res3$FC.B.A <= 1/1.4 & res3$p.B.A < 0.05),])))
  tmp1[2,1] <- sum(grepl("ECOLI",rownames(res3[res3$FC.C.A >= 1.4 & res3$p.C.A < 0.05,])))
  tmp1[2,2] <- sum(grepl("HUMAN",rownames(res3[res3$FC.C.A >= 1.4 & res3$p.C.A < 0.05|(res3$FC.C.A <= 1/1.4 & res3$p.C.A < 0.05),])))
  tmp1[3,1] <- sum(grepl("ECOLI",rownames(res3[res3$FC.D.A >= 1.4 & res3$p.D.A < 0.05,])))
  tmp1[3,2] <- sum(grepl("HUMAN",rownames(res3[res3$FC.D.A >= 1.4 & res3$p.D.A < 0.05|(res3$FC.D.A <= 1/1.4 & res3$p.D.A < 0.05),])))
  tmp1[4,1] <- sum(grepl("ECOLI",rownames(res3[res3$FC.E.A >= 1.4 & res3$p.E.A < 0.05,])))
  tmp1[4,2] <- sum(grepl("HUMAN",rownames(res3[res3$FC.E.A >= 1.4 & res3$p.E.A < 0.05|(res3$FC.E.A <= 1/1.4 & res3$p.E.A < 0.05),])))
  
  tmp1[1,3] <- sum(grepl("ECOLI",rownames(res3[res3$FC.B.A <= 1/1.4 & res3$p.B.A < 0.05,])))
  tmp1[2,3] <- sum(grepl("ECOLI",rownames(res3[res3$FC.C.A <= 1/1.4 & res3$p.C.A < 0.05,])))
  tmp1[3,3] <- sum(grepl("ECOLI",rownames(res3[res3$FC.D.A <= 1/1.4 & res3$p.D.A < 0.05,])))
  tmp1[4,3] <- sum(grepl("ECOLI",rownames(res3[res3$FC.E.A <= 1/1.4 & res3$p.E.A < 0.05,])))
  
  tmp1[,4] <- tmp1[,2]/(tmp1[,1]+tmp1[,2]+tmp1[,3])
  
  rownames(tmp1) <- c("B.A","C.A","D.A","E.A")
  colnames(tmp1) <- c("TP","FP","rev","FADR")
  
  return(tmp1)
}

## add missing value to 
## Function: addMiss (add missing values at desired rate and MNAR ratio)
## MV.rate = overall rate of missing value
## MNAR.ratio = ratio of MNAR (missing not at random) in all missing values
## ini.seed: an initial seed for reproducibility
## version 1.1, bug fixed!(3/29/17)

addMiss <- function(x, MV.rate = 0.3, MNAR.ratio = 0.5, ini.seed = 233){
  
  set.seed(ini.seed)
  seeds <- sample(1:10000,20)
  
  ## Add MNAR missing values according to protein abundance
  
  t_q <- quantile(as.matrix(x), probs = MV.rate)  ## a th quantile of data matrix x
  
  set.seed(seeds[20])
  t_mat <- matrix(rnorm(n = nrow(x)*ncol(x), m = t_q, sd = 0.3),nrow = nrow(x),ncol = ncol(x))  ## threshold matrix
  
  set.seed(seeds[19])
  t_draw <- matrix(rbinom(nrow(x)*ncol(x), 1, MNAR.ratio), nrow = nrow(x),ncol = ncol(x))  ## Bernoulli draw matrix
  
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      if (x[i,j] < t_mat[i,j] & t_draw[i,j] == 1) {
        x[i,j] <- NA
      }
    }
  }  ## add MNAR values
  
  ## sum(is.na(data1_1)) ## number of MNAR missing values
  
  if(sum(is.na(x)) > nrow(x)*ncol(x)*MV.rate) {return(x)}
  else{
    
    ## Add the rest missing values randomly
    
    x2 <- as.vector(as.matrix(x))
    ind1 <- which(is.na(x))
    ind2 <- seq(1:length(x2))
    ind2 <- ind2[!(ind2 %in% ind1)]
    set.seed(seeds[18])
    #ind <- which(x2 %in% sample(x2[!is.na(x2)], as.integer((nrow(x)*ncol(x)*MV.rate-sum(is.na(x))))))
    ind2 <- sample(ind2,as.integer((nrow(x)*ncol(x)*MV.rate-sum(is.na(x)))))
    ind <- sort(c(ind1,ind2))
    x2[ind] <- NA  ##  Add a*b MCAR values randomly
    ## sum(is.na(x2))  ## number of total missing values
    x2 <- matrix(x2, nrow = nrow(x),ncol = ncol(x))
    x2 <- as.data.frame(x2)
    colnames(x2) <- colnames(x)
    rownames(x2) <- rownames(x)
    return(x2)
  }
}

## Function: normImp (imputation using random values from a new normal distribution)
## width: parameter to define sd of new normal distribution = width*(sd of original data)
## down.shift: define mean of new normal distribution = (mean of original data)-down.shift*(sd of original data)
## thres = for protein with real value rate <= thres, NAs will be replaced, set to 1 if all NAs are to be replaced
## group: a factor defines the group information
## ori.seed: set a seed for reproducibility
## Bug fixed 2019/1/23: fixed ave - down.shift*stdev as mean of normal distribution

normImp <- function(x,width = 0.3,down.shift = 2,group = factor(),thres = 1,ori.seed = 666){
  
  ave = mean(as.matrix(x),na.rm = T)
  stdev = sd(as.matrix(x),na.rm = T)
  na.count = sum(is.na(x))
  set.seed(ori.seed)
  res = rnorm(n = na.count,m = ave-down.shift*stdev,sd = stdev*width)
  
  groupSplit <- function(x, group = group) {
    
    df.group <- as.data.frame(table(group))
    df.group$sep <- 1
    
    for(i in 2:dim(df.group)[1]){
      df.group$sep[i] <- df.group$sep[i-1] + df.group$Freq[i-1]
    }
    
    ls <- list()
    for(i in 1:nrow(df.group)) {
      subgroup <- x[,df.group$sep[i]:(df.group$sep[i]+df.group$Freq[i]-1)]
      ls[[i]] <- subgroup
    }
    
    names(ls) <- df.group$Var1
    return(ls)
  }
  
  lib <- groupSplit(x,group = group)
  idx = matrix(F,nrow = nrow(x),ncol = ncol(x))  ## determine which NA needs to be replaced
  lib_idx = groupSplit(as.data.frame(idx),group = group)
  
  ## assign TRUE to NAs to be replaced according to threshold
  for(i in 1:length(lib)){
    tmp = lib[[i]]
    tmp_idx = lib_idx[[i]]
    t = ncol(tmp)*(1-thres)
    for(j in 1:nrow(tmp)){
      if(sum(is.na(tmp[j,])) >= t){
        tmp_idx[j,] <- is.na(tmp[j,])
      }
    }
    lib_idx[[i]] <- tmp_idx
  }
  res_idx <- lib_idx[[1]]  ## final index matrix
  for(i in 2:(length(lib_idx))){res_idx <- cbind(res_idx,lib_idx[[i]])}
  
  ## replace NAs based on index matrix with random values from random distribution
  c=1
  for(i in 1:nrow(x)){
    for(j in 1:ncol(x)){
      if(res_idx[i,j]){x[i,j]=res[c];c=c+1}
    }
  }
  
  return(x)
}

## function: groupSplit (split a dataframe by factor and by column, return a list)
## group = factor
## v1.1 bug fixed

groupSplit <- function(x, group = group) {
  
  df.group <- as.data.frame(table(group))
  df.group$sep <- 1
  
  for(i in 2:dim(df.group)[1]){
    df.group$sep[i] <- df.group$sep[i-1] + df.group$Freq[i-1]
  }
  
  ls <- list()
  for(i in 1:nrow(df.group)) {
    subgroup <- x[,df.group$sep[i]:(df.group$sep[i]+df.group$Freq[i]-1)]
    ls[[i]] <- subgroup
  }
  
  names(ls) <- df.group$Var1
  return(ls)
}

## Function: naCon (return number of conditions with no observed values)
## group = factor with condition information

naCon <- function(x, group = group) {
  
  c <- 0
  df.group <- as.data.frame(table(group))
  df.group$sep <- 1
  
  for(i in 2:nrow(df.group)){
    df.group$sep[i] <- df.group$sep[i-1] + df.group$Freq[i]
  }
  
  for(i in 1:nrow(df.group)){
    
    if(i+1 <= nrow(df.group)){
      if(sum(is.na(x[df.group$sep[i]:((df.group$sep[i+1])-1)])) == df.group$Freq[i]){c <- c+1}
    }
    if(i+1 > nrow(df.group)){
      if(sum(is.na(x[df.group$sep[i]:length(x)])) == df.group$Freq[i]){c <- c+1}
    }
  }
  return(c)
}

## function:lloq3  (incorporate function:groupSplit)
## This function also returns T/F to index imputed value
## method = c("global","conditional"), missing value rate will be calculated accordingly
## thres: for protein with real value rate <= thres, NAs will be replaced by lloq
## adj = adjustment parameter that will be used to substract from min.value for lloq replacement (added 3/13/17)

lloq3 <- function(x, group = factor(), thres = 0.4, adj = 0){
  #  if(!method %in% c("global","conditional")){stop("Wrong method! Please choose from global and conditional.")}
  
  min.value = min(abs(x), na.rm = T)-adj
  
  groupSplit <- function(x, group = group) {
    
    df.group <- as.data.frame(table(group))
    df.group$sep <- 1
    
    for(i in 2:dim(df.group)[1]){
      df.group$sep[i] <- df.group$sep[i-1] + df.group$Freq[i-1]
    }
    
    ls <- list()
    for(i in 1:nrow(df.group)) {
      subgroup <- x[,df.group$sep[i]:(df.group$sep[i]+df.group$Freq[i]-1)]
      ls[[i]] <- subgroup
    }
    
    names(ls) <- df.group$Var1
    return(ls)
  }
  
  lib = groupSplit(x,group = group)
  idx = matrix(F,nrow = nrow(x),ncol = ncol(x))
  #colnames(idx) = paste0("imp_",colnames(x))
  lib_idx = groupSplit(as.data.frame(idx),group = group)
  
  for(i in 1:length(lib)){
    tmp = lib[[i]]
    tmp_idx = lib_idx[[i]]
    t = ncol(tmp)*(1-thres)
    for(j in 1:nrow(tmp)){
      if(sum(is.na(tmp[j,])) >= t){
        tmp_idx[j,] <- is.na(tmp[j,])
        v = tmp[j,]
        v[is.na(v)] <- min.value
        tmp[j,] <- v
      }
    }
    lib[[i]] <- tmp
    lib_idx[[i]] <- tmp_idx
  }
  
  res <- lib[[1]]
  for(i in 2:(length(lib))){res <- cbind(res,lib[[i]])}
  res_idx <- lib_idx[[1]]
  for(i in 2:(length(lib_idx))){res_idx <- cbind(res_idx,lib_idx[[i]])}
  
  res <- as.data.frame(res)
  colnames(res) <- colnames(x)
  rownames(res) <- rownames(x)
  res_idx <- as.data.frame(res_idx)
  colnames(res_idx) <- paste0("imp_",colnames(x))
  rownames(res_idx) <- rownames(x)
  res2 <- cbind(res,res_idx)
  return(res2)
}

## function:lloq  (requires function:groupSplit)
## method = c("global","conditional"), missing value rate will be calculated accordingly
## thres: for protein with real value rate <= thres, NAs will be replaced by lloq
## adj = adjustment parameter that will be used to substract from min.value for lloq replacement (added 3/13/17)

lloq <- function(x, group = factor(),method = "global", thres = 0.3, adj = 0){
  if(!method %in% c("global","conditional")){stop("Wrong method! Please choose from global and conditional.")}
  
  min.value <- min(abs(x), na.rm = T)-adj
  
  minReplace <- function(y,t = thres, m = min.value){
    if(sum(is.na(y)) > length(y)*(1-t)){y[is.na(y)] <- m}
    return(y)
  }
  
  if(method == "global"){   
    res <- apply(as.matrix(x),1,minReplace)
    res <- t(res)
  }
  if(method == "conditional"){
    lib <- groupSplit(x,group = group)
    for(i in 1:length(lib)){
      lib[[i]] <- apply(as.matrix(lib[[i]]),1,minReplace)
      lib[[i]] <- t(lib[[i]])
    }
    res <- vector()
    for(i in 1:length(lib)){res <- cbind(res,lib[[i]])}
  }
  return(res)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## Function: missCut2(v2.0): requires function groupSplit
## missing.type = c("zero","NA")
## cut.g = global threshold, protein ID with real value rate over cut.g will be kept
## cut.c = conditional threshold, protein ID with real value rate over cut.c in >= one condition will be kept

missCut2 <- function(x,
                     group = factor(),
                     missing.type = "NA",
                     global.cut = T,
                     condition.cut = T,
                     cut.g = 0.5,
                     cut.c = 0.6){
  
  if(missing.type == "zero"){x[x == 0] <- NA}
  c_glo <- vector()
  c_con <- vector()
  if(global.cut){
    x_g <- x[rowSums(is.na(x)) <= ncol(x)*(1-cut.g),]
    c_glo <- rownames(x_g)
  }
  if(condition.cut){
    lib <- groupSplit(x,group = group)
    for(i in 1:length(lib)){
      tmp <- lib[[i]][rowSums(is.na(lib[[i]])) <= ncol(lib[[i]])*(1-cut.c),]
      c_con <- c(c_con,rownames(tmp))
    }
  }
  pro_keep <- unique(c(c_glo,c_con))
  x_keep <- x[pro_keep,]
  p1 <- paste("In total", nrow(x),"proteins,", nrow(x_keep), "are kept after missing value cut off.")
  print(p1)
  if(missing.type == "zero"){x_keep[is.na(x_keep)] <- 0}
  return(x_keep)
}

## function: groupSplit (split a dataframe by factor and by column, return a list)
## group = factor
## v1.1 bug fixed

groupSplit <- function(x, group = group) {
  
  df.group <- as.data.frame(table(group))
  df.group$sep <- 1
  
  for(i in 2:dim(df.group)[1]){
    df.group$sep[i] <- df.group$sep[i-1] + df.group$Freq[i-1]
  }
  
  ls <- list()
  for(i in 1:nrow(df.group)) {
    subgroup <- x[,df.group$sep[i]:(df.group$sep[i]+df.group$Freq[i]-1)]
    ls[[i]] <- subgroup
  }
  
  names(ls) <- df.group$Var1
  return(ls)
}