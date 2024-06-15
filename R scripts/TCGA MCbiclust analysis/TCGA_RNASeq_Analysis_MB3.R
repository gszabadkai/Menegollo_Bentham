# Analysis of the TCGA RNASeq data using special 
# MCbiclust adapted function to find bicluster matching
# regulation previously found MB2 and MB3. 
# TCGA RNASeq data in this analysis has had the normal 
# samples removed.


# 1. Load BRCA RNASeq data ---
library(MCbiclust)
library(ggplot2)
library(WriteXLS)
library(dplyr)

# Data from previous analysis loaded in here
load('../METABRIC_starting_data_28_01_2021_corr_groups.Rdata')
cor_vec_all <- mito.cv[[1]]

# Select the top 1000 genes from the BRCA microarray correlation vector
top1000.BRCA <- order(abs(cor_vec_all),decreasing = T)[c(1:1000)]

# Split these genes into two groups based on 
# whether they have positive/negative correlation
gene.group1.loc <- top1000.BRCA[which(cor_vec_all[top1000.BRCA] > 0)]
gene.group2.loc <- top1000.BRCA[which(cor_vec_all[top1000.BRCA] < 0)]

gene.group1 <- METABRIC.data[gene.group1.loc,1]
gene.group2 <- METABRIC.data[gene.group2.loc,1]

# Load RNASeq data
load("Legion_BRCA_RNAseq1.RData")
#RNA_seq_v2_mat2 <- log2(RNA_seq_v2_mat + 1)
# Find the corresponding gene groups in the RNAseq dataset

RNAseq.group1 <- which(row.names(RNA_seq_v2_mat) %in% gene.group1)
RNAseq.group2 <- which(row.names(RNA_seq_v2_mat) %in% gene.group2)

# Set random seed to make work reproducible
set.seed(1)

# Using the alternative bicluster find the seed. ----
FindSeedGroups <-function(gem, seed.size, iterations, group1.loc, group2.loc,
                          initial.seed = NULL, full.detail = FALSE){
  gem <- gem[c(group1.loc,group2.loc),]
  sample.size <- dim(gem)[2]
  
  if (length(initial.seed) == seed.size){
    main.subsamp <- initial.seed}
  else {
    main.subsamp <- sample(seq(length = sample.size), seed.size)}
  
  gem.t<-t(gem)
  
  # If rows contain only zero/constant values correlation results in NA values
  zero.row.update <- FSGZeroRowUpdate(gem,main.subsamp,group1.loc,group2.loc)
  
  test.cor.mat <- zero.row.update[[1]]
  group1.loc.upd <- zero.row.update[[2]]
  group2.loc.upd <- zero.row.update[[3]]
  
  test.cor.mat1 <- test.cor.mat[group1.loc.upd, group1.loc.upd]
  test.cor.mat2 <- test.cor.mat[group2.loc.upd, group2.loc.upd]
  test.cor.mat3 <- test.cor.mat[group1.loc.upd, group2.loc.upd]
  
  # Normalise to group length (corresponding to area on heatmap)
  ta1 <- mean(test.cor.mat1) / (length(group1.loc) * length(group1.loc))
  ta2 <- mean(test.cor.mat2) / (length(group2.loc) * length(group2.loc))
  ta3 <- (-mean(test.cor.mat3) / (length(group1.loc) * length(group2.loc)))
  
  test.val <- mean(c(ta1,ta2,ta3))
  remove.sample <- sample(seq(length = seed.size), iterations, replace=TRUE)
  
  if(full.detail == TRUE){
    sample.list <- list()
    sample.list[[1]] <- c(0, test.val, mean(abs(test.cor.mat)), main.subsamp)
    sample.rank <- 1
  }
  
  for(i in 1:iterations){
    # remove one sample randomly
    test.replace <- sample(seq(length = sample.size)[-c(main.subsamp)],1)
    test.subsamp <- main.subsamp[-remove.sample[i]]
    # replace with a new sample
    test.subsamp <- c(test.subsamp, test.replace)
    
    # retest zero.rows
    zero.row.update <- FSGZeroRowUpdate(gem,test.subsamp,group1.loc,group2.loc)
    
    test.cor.mat <- zero.row.update[[1]]
    group1.loc.upd <- zero.row.update[[2]]
    group2.loc.upd <- zero.row.update[[3]]
    
    test.cor.mat1 <- test.cor.mat[group1.loc.upd,group1.loc.upd]
    test.cor.mat2 <- test.cor.mat[group2.loc.upd,group2.loc.upd]
    test.cor.mat3 <- test.cor.mat[group1.loc.upd,group2.loc.upd]
    
    # Normalise to group length (corresponding to area on heatmap)      
    ta1 <- mean(test.cor.mat1) / (length(group1.loc) * length(group1.loc))     
    ta2 <- mean(test.cor.mat2) / (length(group2.loc) * length(group2.loc))      
    ta3 <- (-mean(test.cor.mat3) / (length(group1.loc) * length(group2.loc)))
    
    test.val2 <- mean(c(ta1, ta2, ta3))
    
    if(test.val2 > test.val){
      main.subsamp <- test.subsamp
      test.val <- test.val2
      print.val <- mean(abs(test.cor.mat))
      if(full.detail == TRUE){
        sample.rank <- sample.rank + 1
        sample.list[[sample.rank]] <- c(i, test.val, print.val, main.subsamp)
      }
    }
    if(i %% 100 == 0){
      print(c(i,test.val,print.val))}
  }
  return(main.subsamp)
}


FSGZeroRowUpdate <- function(gem,seed,group1.loc,group2.loc){
  gem.t <- t(gem)
  zero.rows <- which(apply(X = gem[, seed],MARGIN = 1,FUN = sd) == 0)
  
  if(length(zero.rows) != 0){
    test.cor.mat <- cor(gem.t[seed, -zero.rows])
    
    group1.zr <- which(zero.rows <= length(group1.loc))
    group2.zr <- which(zero.rows > length(group1.loc))
    
    if(length(group1.zr) > 0){
      g1 <- length(group1.loc) - length(group1.zr)
      group1.loc.upd <- group1.loc[1:g1]}
    else{
      group1.loc.upd <- group1.loc}
    
    if(length(group2.zr) > 0){
      g2 <- length(group2.loc) - length(group2.zr)
      group2.loc.upd <- group2.loc[1:g2] - length(group1.zr)}
    else{
      group2.loc.upd <- group2.loc - length(group1.zr)}
  }
  else{
    test.cor.mat <- cor(t(gem[ ,seed]))
    group1.loc.upd <- group1.loc
    group2.loc.upd <- group2.loc}
  return(list(test.cor.mat,group1.loc.upd,group2.loc.upd))
}

# Remove normal samples ----

sample.codes <- substr(colnames(RNA_seq_v2_mat),14,15)
# remove normal (11) and metastatic (06)
s1 <- which(sample.codes %in% c('11','06'))

RNA.ma.pat <- log2(RNA_seq_v2_mat[c(RNAseq.group1,RNAseq.group2), -s1] + 1)
RNA.ma.pat2 <- RNA.ma.pat - rowMeans(RNA.ma.pat)

ma.len1 <- seq(length = length(RNAseq.group1))
ma.len2 <- seq(length = dim(RNA.ma.pat)[1])[-seq(length = length(RNAseq.group1))]

RNA.seed.groups <- FindSeedGroups(gem = RNA.ma.pat,
                                  seed.size = 10,
                                  iterations = 10000,
                                  group1.loc = ma.len1,
                                  group2.loc = ma.len2)

zero.rows <- which(rowMeans(RNA.ma.pat[,RNA.seed.groups]) == 0)

# Sort samples ----
new.sort <- SampleSort(gem = RNA.ma.pat[-zero.rows,], seed = RNA.seed.groups,
                       num.cores = 1)
  
cv.test <- CVEval(RNA.ma.pat[-zero.rows,], log2(RNA_seq_v2_mat[, -s1] + 1),
                  seed = new.sort[seq(10)],splits = 8)

# Find TCGA pattern + fork ----

new.pc1 <- PC1VecFun(RNA.ma.pat2[-zero.rows,],new.sort,n = 10)

# Need to make df with PAM50 info ----
load("PAM50results.RData")
RNAseq.pam50 <- read.delim("RNAseq_log_pam50scores.txt")


# note that labels were loaded wrong from PAM50 - missing 1st and no label on last
RNAseq.pam50$X <- colnames(RNA_seq_v2_mat)

new.MB3.sort.df <- data.frame(X = colnames(RNA.ma.pat)[new.sort],
                          MB3.pc1 = new.pc1, MB3.index =  seq_len(length(new.sort)),
                          PAM50 = (RNAseq.pam50$Call)[-s1][new.sort],
                          LumA = (RNAseq.pam50$LumA)[-s1][new.sort],
                          LumB = (RNAseq.pam50$LumB)[-s1][new.sort])


# WriteXLS(new.sort.df, 'TCGA_RNASeq_nonorm_MB3.xlsx')

TCGA.MB3.cv <- CVEval(RNA.ma.pat[-zero.rows,], log2(RNA_seq_v2_mat[,-s1] + 1),
                  new.sort[seq(10)],splits = 8)

TCGA.MB3.CV.RNAseq.df <- data.frame(Genes = rownames(RNA_seq_v2_mat),
                                CV = as.numeric(TCGA.MB3.cv))


TCGA.MB3.RNAseq.df <- new.MB3.sort.df 
save(list = c('TCGA.MB3.RNAseq.df','TCGA.MB3.CV.RNAseq.df'),
     file = 'TCGA_MB3_RNASeq_data_nonorm.RData')



# MB3 lookup ----
ggplot(new.MB3.sort.df , aes(MB3.index, MB3.pc1)) +
  geom_point(aes(colour = PAM50))
ggplot(new.MB3.sort.df , aes(MB3.index, MB3.pc1)) +
  geom_point(aes(colour = LumA))
ggplot(new.MB3.sort.df , aes(MB3.index, MB3.pc1)) +
  geom_point(aes(colour = LumB))
ggplot(new.MB3.sort.df , aes(MB3.index, MB3.pc1)) +
  geom_point(aes(colour = LumA-LumB))

