# compile all data and derived parameters from METABRIC data
# !!!! new nomenclature massive bicluster MB1 = ICT1; MB2 = Mito2; MB3 = Mito1 




library(dplyr)
library(gplots)
library(tidyr)
library(tibble)
library(stringr)
library(MCbiclust)
library(gdata)
library(WriteXLS)
library(biomaRt)
library(propagate)
library(GSVA)
library(readxl)


# Load data----
load('METABRIC_DATA.RData')
load('METABRIC_gene_lists.RData')
load('METABRIC_sort_data.RData')
load('METABRIC_PC1_GSEA.RData')

# Make dataframe for fork parameters pc1 and index, fork and forkscale ----
sort.df <- data.frame(sample = Complete_METABRIC_Clinical_Features_Data[,1],
                      MB3.index = order(mito.sort[[1]]),
                      MB3.pc1  = mito.pc1[[1]][order(mito.sort[[1]])],
                      MB2.index = order(mito.sort[[2]]),
                      MB2.pc1  = mito.pc1[[2]][order(mito.sort[[2]])],
                      Mito3.index = order(mito.sort[[3]]),
                      Mito3.pc1  = mito.pc1[[3]][order(mito.sort[[3]])],
                      MB1.index = order(ICT1.sort[[1]]),
                      MB1.pc1  = ICT1.pc1[[1]][order(ICT1.sort[[1]])],
                      Random1.index = order(random.sort[[1]]),
                      Random1.pc1  = random.pc1[[1]][order(random.sort[[1]])],
                      Random2.index = order(random.sort[[2]]),
                      Random2.pc1  = random.pc1[[2]][order(random.sort[[2]])])


all.clinical.df <- left_join(sort.df,
                             Complete_METABRIC_Clinical_Features_Data,
                             by = 'sample')

colnames(all.clinical.df)[35] <- 'PAM50'
colnames(all.clinical.df)[19] <- 'ER status'
colnames(all.clinical.df)[18] <- "Hist_type"

all.clinical.df <- mutate(all.clinical.df, `Hist type` = ifelse(grepl("IDC", `Hist_type`), "IDC",
                                                         ifelse(grepl("ILC", `Hist_type`), "ILC",
                                                        ifelse(grepl("DCIS", `Hist_type`), "DCIS", "Other"))))


all.clinical.df <- mutate(all.clinical.df, MB1.forkscale = MB1.pc1/MB1.index)
all.clinical.df <- mutate(all.clinical.df, MB3.pc1.rev = MB3.pc1 * -1)
all.clinical.df <- mutate(all.clinical.df, MB3.forkscale = MB3.pc1.rev/MB3.index)
all.clinical.df <- mutate(all.clinical.df, MB2.forkscale = MB2.pc1/MB2.index)

all.clinical.df <- mutate(all.clinical.df, MB1.forkscale.fork_0.2 = ifelse(`MB1.forkscale` <= -0.02 , "Lower", ifelse(`MB1.forkscale` >= 0.02 , "Upper", "none")))
all.clinical.df <- mutate(all.clinical.df, MB1.forkscale.fork_0.4 = ifelse(`MB1.forkscale` <= -0.04 , "Lower", ifelse(`MB1.forkscale` >= 0.04 , "Upper", "none")))
all.clinical.df <- mutate(all.clinical.df, MB1.forkscale.fork_0.8 = ifelse(`MB1.forkscale` <= -0.08 , "Lower", ifelse(`MB1.forkscale` >= 0.08 , "Upper", "none")))
all.clinical.df <- mutate(all.clinical.df, MB2.forkscale.fork_0.2 = ifelse(`MB2.forkscale` <= -0.02 , "Lower", ifelse(`MB2.forkscale` >= 0.02 , "Upper", "none")))
all.clinical.df <- mutate(all.clinical.df, MB2.forkscale.fork_0.4 = ifelse(`MB2.forkscale` <= -0.04 , "Lower", ifelse(`MB2.forkscale` >= 0.04 , "Upper", "none")))
all.clinical.df <- mutate(all.clinical.df, MB2.forkscale.fork_0.8 = ifelse(`MB2.forkscale` <= -0.08 , "Lower", ifelse(`MB2.forkscale` >= 0.08 , "Upper", "none")))
all.clinical.df <- mutate(all.clinical.df, MB3.forkscale.fork_0.2 = ifelse(`MB3.forkscale` <= -0.02 , "Lower", ifelse(`MB3.forkscale` >= 0.02 , "Upper", "none")))
all.clinical.df <- mutate(all.clinical.df, MB3.forkscale.fork_0.4 = ifelse(`MB3.forkscale` <= -0.04 , "Lower", ifelse(`MB3.forkscale` >= 0.04 , "Upper", "none")))
all.clinical.df <- mutate(all.clinical.df, MB3.forkscale.fork_0.8 = ifelse(`MB3.forkscale` <= -0.08 , "Lower", ifelse(`MB3.forkscale` >= 0.08 , "Upper", "none")))


all.clinical.df$MB1.forkscale.fork_0.2 <- paste("MB1",all.clinical.df$MB1.forkscale.fork_0.2,sep="_")
all.clinical.df$MB2.forkscale.fork_0.2 <- paste("MB2",all.clinical.df$MB2.forkscale.fork_0.2,sep="_")
all.clinical.df$MB3.forkscale.fork_0.2 <- paste("MB3",all.clinical.df$MB3.forkscale.fork_0.2,sep="_")
all.clinical.df$MB1.forkscale.fork_0.4 <- paste("MB1",all.clinical.df$MB1.forkscale.fork_0.4,sep="_")
all.clinical.df$MB2.forkscale.fork_0.4 <- paste("MB2",all.clinical.df$MB2.forkscale.fork_0.4,sep="_")
all.clinical.df$MB3.forkscale.fork_0.4 <- paste("MB3",all.clinical.df$MB3.forkscale.fork_0.4,sep="_")
all.clinical.df$MB1.forkscale.fork_0.8 <- paste("MB1",all.clinical.df$MB1.forkscale.fork_0.8,sep="_")
all.clinical.df$MB2.forkscale.fork_0.8 <- paste("MB2",all.clinical.df$MB2.forkscale.fork_0.8,sep="_")
all.clinical.df$MB3.forkscale.fork_0.8 <- paste("MB3",all.clinical.df$MB3.forkscale.fork_0.8,sep="_")

# plot(density(all.clinical.df2$Mito1.forkscale.log))
# min(all.clinical.df2$ICT1.forkscale.log)


# fork classifications based on ThresholdBic: Methods for defining a bicluster (https://rdrr.io/bioc/MCbiclust/man/ThresholdBic.html) ('Bobby's method') -----

MB1.thr <- ThresholdBic(ICT1.cv[[1]], sort.order = ICT1.sort[[1]],
                         pc1 = ICT1.pc1[[1]],samp.sig = 0.8)

MB1.fork <- ForkClassifier(mito.pc1[[1]],
                            samp.num = length(MB1.thr[[2]]))
all.clinical.df$MB1.Thresholdfork <- MB1.fork[order(mito.sort[[1]])]


mito.pc1[[4]] <- mito.pc1[[1]] * -1

MB3.thr <- ThresholdBic(mito.cv[[1]], sort.order = mito.sort[[1]],
                          pc1 = mito.pc1[[4]],samp.sig = 0.95)
MB3.fork <- ForkClassifier(mito.pc1[[4]],
                             samp.num = length(MB3.thr[[2]]))
all.clinical.df$MB3.Thresholdfork <- MB3.fork[order(mito.sort[[1]])]


MB2.thr <- ThresholdBic(mito.cv[[2]], sort.order = mito.sort[[2]],
                          pc1 = mito.pc1[[2]],samp.sig = 0.95)
MB2.fork <- ForkClassifier(mito.pc1[[2]],
                             samp.num = length(MB2.thr[[2]]))
all.clinical.df$MB2.Thresholdfork <- MB2.fork[order(mito.sort[[2]])]

all.clinical.df$MB1.Thresholdfork <- paste('MB1',all.clinical.df$MB1.Thresholdfork,sep = '_')
all.clinical.df$MB3.Thresholdfork <- paste('MB3',all.clinical.df$MB3.Thresholdfork,sep = '_')
all.clinical.df$MB2.Thresholdfork <- paste('MB2',all.clinical.df$MB2.Thresholdfork,sep = '_')
all.clinical.df[,c(44:55)] <- lapply(all.clinical.df[,c(44:55)], function(x) {gsub("Lower", "LF", x)})
all.clinical.df[,c(44:55)] <- lapply(all.clinical.df[,c(44:55)], function(x) {gsub("Upper", "UF", x)})
all.clinical.df[,c(53:55)] <- lapply(all.clinical.df[,c(53:55)], function(x) {gsub("None", "none", x)})


# the best forkscale log: pc1/log(index) created for all three biclusters -----
all.clinical.df <- all.clinical.df %>% 
  mutate(MB1.forkscale.log = MB1.pc1/(log(MB1.index))) %>% 
  mutate(MB2.forkscale.log = MB2.pc1/(log(MB2.index))) %>% 
  mutate(MB3.forkscale.log = -(MB3.pc1/(log(MB3.index))))



# Save METABRIC_data_w_forkscales.xls ----
save(list = c('all.clinical.df'),file = 'METABRIC_data_w_forkscales.xls')


# Add ESTIMATE results----
# ## Using the ESTIMATE package to calculate immune and stromal score for the METABRIC data
# 
# See <http://bioinformatics.mdanderson.org/main/ESTIMATE:Overview> for details of the estimate package.
# 
# ESTIMATE (Estimation of STromal and Immune cells in MAlignant Tumor tissues using Expression data) is a tool for predicting tumour purity, and the presence of infiltrating stromal/immune cells in tumour tissues using gene expression data. ESTIMATE algorithm is based on single sample Gene Set Enrichment Analysis and generates three scores:
#   
# 1. Stromal score (that captures the presence of stroma in tumour tissue)
# 2. Immune score (that represents the infiltration of immune cells in tumour tissue), and
# 3. estimate score (that infers tumour purity)


#  original location: tmp.df <- read.delim("~/Google Drive/My Drive/data/CCLE and BRCA project from 2013/METABRIC/R files/METABRIC_estimate_score.gct", comment.char="#", skip = 1)

tmp.df <- read.delim("METABRIC_estimate_score.gct", comment.char="#", skip = 1)


METABRIC.estimate <- data.frame(
  sample = as.character(sapply(tmp.df[1,-c(1,2)], as.character)),
  StromalScore = as.numeric(sapply(tmp.df [2,-c(1,2)],FUN = as.character)),
  ImmuneScore = as.numeric(sapply(tmp.df [3,-c(1,2)],FUN = as.character)),
  ESTIMATEScore = as.numeric(sapply(tmp.df [4,-c(1,2)],FUN = as.character))
)

METABRIC.estimate$sample <- gsub('\\.','-',METABRIC.estimate$sample)

all.clinical.df2 <- left_join(all.clinical.df, METABRIC.estimate,"sample")

# Add copy number data ---- not used in paper

# # METABRIC Copy number
# 
# # Load Data
# #  original location:  load('~/Google Drive/My Drive/R/METABRIC copy number/Complete_METABRIC_Copy_Number_Data.rbin')
# 
# load('Complete_METABRIC_Copy_Number_Data.rbin')
# 
# METABRIC.copy <- as.data.frame(Complete_METABRIC_Copy_Number_Data@assayData$exprs)
# 
# # Get gene names (eg = entrez gene?)
# library(biomaRt)
# 
# # if (!requireNamespace("BiocManager", quietly = TRUE))
# #   install.packages("BiocManager")
# # BiocManager::install("biomaRt", version = "3.8")
# 
# ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
# attributes = listAttributes(ensembl)
# 
# probe2gene <- getBM(attributes=c('entrezgene_id', 'external_gene_name',
#                                  'chromosome_name','start_position','band'),
#                     filters = 'entrezgene_id',
#                     values = gsub('_eg','',row.names(METABRIC.copy)),
#                     mart = ensembl)
# 
# chr.loc <- c(grep('CHR',unique(probe2gene$chromosome_name)),
#              grep('\\.',unique(probe2gene$chromosome_name)))
# 
# probe2gene2 <- probe2gene %>%
#   filter(!chromosome_name %in% unique(probe2gene$chromosome_name)[chr.loc]) %>%
#   mutate(cytoband = paste(chromosome_name, band,sep = ''))
# 
# probe2gene2 <- arrange(probe2gene2, chromosome_name, start_position)
# 
# new.order.vec <- c(as.character(c(2:22)), 'X', 'Y')
# new.order <- which(probe2gene2$chromosome_name == '1')
# for(i in 1:23){
#   new.order <- c(new.order,which(probe2gene2$chromosome_name == new.order.vec[i]))
# }
# 
# probe2gene2 <- probe2gene2[new.order,]
# 
# probe2gene$entrezgene <- as.character(probe2gene$entrezgene)
# # Remove entrez genes with multiple gene names/starting positions
# unique.loc <- which(sapply(seq(18026),
#                            FUN = function(x) length(which(probe2gene2$entrezgene ==
#                                                             (unique(probe2gene2$entrezgene))[x]))) > 1)
# 
# probe2gene2$entrezgene <- as.character(probe2gene2$entrezgene)
# probe2gene2 <- filter(probe2gene2, !(entrezgene %in% unique(probe2gene$entrezgene)[unique.loc]))
# 
# 
# METABRIC.copy$entrezgene <- gsub('_eg','',row.names(METABRIC.copy))
# 
# METABRIC.copy2 <- inner_join(probe2gene2, METABRIC.copy, by = "entrezgene")


# 
# Getting gene sets----

# the *.cv files contain cv values for all genes accross the samples, showing Pearson's correlation with the genes in the bicluster
# ICT1.hicor.loc and mito.genes.loc vectors of genes used as the geneset to start MCbiclust. 
# ICT1.hicor.loc is generated by selecting 1000 genes showing highest correlation with the ICT1 gene (MRPL59), mito.genes.loc are the genes in Mitocarta 1.0 
# dataframes are created with the names and cv values of these genes. Two groups are also created by clustering the potivie and negatively correlating genes
# calculation of ICT1.hicor.loc is based on the correlation of genes in the dataset with ICT1=MRPL58 gene, example below is for another dataset (see NAR paper), below we use ICT1.hicor.loc determined previously by Bobby for the METABRIC dataset
# ICT1.loc <- which(row.names(combat_edata2) == "MRPL58")
# ICT1.cv <- sapply(seq(dim(combat_edata2)[1]),
#                   FUN= function(x) cor(combat_edata2[ICT1.loc,],combat_edata2[x,]))
# ICT1.hicor.loc <- order(abs(ICT1.cv),decreasing = T)[seq(1000)]

###MB1 (former ICT1) ----
# 1. get genes and their CVs in the bicluster:

# MB1.bicluster.genes <- data.frame(Gene=(METABRIC.data[,1][MB1.hicor.loc]), cv=as.vector(MB1.cv[[1]])[MB1.hicor.loc], stringsAsFactors = FALSE) 

MB1.bicluster.genes <- data.frame(Gene=(METABRIC.data[ICT1.hicor.loc,1]), cv=as.vector(ICT1.cv[[1]])[ICT1.hicor.loc], stringsAsFactors = FALSE) 

# 2. get genes with the highest CV from the whole genome:
MB1.hi.cv.loc <- order(abs(ICT1.cv[[1]]),decreasing = TRUE)[seq(1000)]
MB1.hi.cv.genes <- data.frame(Gene = (METABRIC.data[,1][MB1.hi.cv.loc]), cv=as.vector(ICT1.cv[[1]])[MB1.hi.cv.loc], stringsAsFactors = FALSE) 

# 3. find the two anti-correlating groups by clustering the genes in the top 10 samples by heatmap.2
# a. for genes with the highest CV from the whole genome:

# pdf("MB1.all.cv.heat.pdf")

MB1.all.cv.heat <- heatmap.2(cor(t(METABRIC.data[,-1][MB1.hi.cv.loc,ICT1.sort[[1]][seq(10)]])),
                             trace = 'none')



# dev.off()


# b. for the highest correlating genes in the bicluster (note that here instead of simple cor we use HclustGenesHiCor which trims the dendrogram to the highest correlating groups)

MB1.hicor.genes <- as.numeric(HclustGenesHiCor(METABRIC.data[,-1][ICT1.hicor.loc,],
                                                 ICT1.top.seed[[1]],
                                                 cuts = 8))

# pdf("MB1.biclustergenes.cv.heat.pdf")

MB1.biclustergenes.cv.heat <- heatmap.2(cor(t(METABRIC.data[,-1][ICT1.hicor.loc[MB1.hicor.genes],
                                                           ICT1.sort[[1]][seq(10)]])),
                                  trace = 'none')
# dev.off()

# metabolic genes

metab.genes.xls <- read_excel('GS metabolic genes list.xlsx', sheet=1)
metab.loc <- which(METABRIC.data[,1] %in% as.data.frame(metab.genes.xls)[,2])
metab.genes <- data.frame(Gene=(METABRIC.data[,1][metab.loc]), cv=as.vector(ICT1.cv[[1]])[metab.loc], stringsAsFactors = FALSE)
  # metabolic genes highly correlating with the bicluster
MB1.hicor.metab.genes <- as.numeric(HclustGenesHiCor(METABRIC.data[,-1][metab.loc,],
                                                ICT1.top.seed[[1]],
                                                cuts = 6))

# pdf("MB1.metab.genes.cv.heat.pdf")

MB1.metab.genes.cv.heat <- heatmap.2(cor(t(METABRIC.data[,-1][metab.loc[MB1.hicor.metab.genes],
                                                                  ICT1.sort[[1]][seq(10)]])),
                                         trace = 'none')

# dev.off()


# get genes together in dataframes with names and CVs in a list from all the above

df.list <- list()

df.list$MB1.hi.cv.genes.cv <- as.data.frame(MB1.hi.cv.genes) # genes with the highest CV from the whole genome:
df.list$MB1.bicluster.genes.cv <-  as.data.frame(MB1.bicluster.genes) # genes and their CVs in the bicluster
df.list$MB1.allgenome.cv <- data.frame(Gene = METABRIC.data[,1],CV = ICT1.cv[[1]]) #CVs of all genes in the genome
df.list$MB1.all.metab.genes.cv <- data.frame(Gene = as.character(metab.genes[,1]), CV = df.list$MB1.allgenome.cv[match(metab.genes[,1], df.list$MB1.allgenome.cv$Gene),2], stringsAsFactors = FALSE) #CVs of all metabolic genes in the genome
df.list$MB1.bicluster.group1 <- data.frame(Gene = METABRIC.data[labels(MB1.biclustergenes.cv.heat$rowDendrogram[[2]]),1], cv = MB1.bicluster.genes[match(c(METABRIC.data[labels(MB1.biclustergenes.cv.heat$rowDendrogram[[2]]),1]),MB1.bicluster.genes$Gene),2], stringsAsFactors = FALSE)
df.list$MB1.bicluster.group2 <- data.frame(Gene = METABRIC.data[labels(MB1.biclustergenes.cv.heat$rowDendrogram[[1]]),1], cv = MB1.bicluster.genes[match(c(METABRIC.data[labels(MB1.biclustergenes.cv.heat$rowDendrogram[[1]]),1]),MB1.bicluster.genes$Gene),2], stringsAsFactors = FALSE)
df.list$MB1.hi.cv.group1 <- data.frame(Gene = METABRIC.data[labels(MB1.all.cv.heat$rowDendrogram[[1]]),1], cv = MB1.hi.cv.genes[match(c(METABRIC.data[labels(MB1.all.cv.heat$rowDendrogram[[1]]),1]),MB1.hi.cv.genes$Gene),2], stringsAsFactors = FALSE)
df.list$MB1.hi.cv.group2 <- data.frame(Gene = METABRIC.data[labels(MB1.all.cv.heat$rowDendrogram[[2]]),1], cv = MB1.hi.cv.genes[match(c(METABRIC.data[labels(MB1.all.cv.heat$rowDendrogram[[2]]),1]),MB1.hi.cv.genes$Gene),2], stringsAsFactors = FALSE)
df.list$MB1.metab.group1 <- data.frame(Gene = METABRIC.data[labels(MB1.metab.genes.cv.heat$rowDendrogram[[2]]),1], cv = metab.genes[match(c(METABRIC.data[labels(MB1.metab.genes.cv.heat$rowDendrogram[[2]]),1]),metab.genes$Gene),2], stringsAsFactors = FALSE)
df.list$MB1.metab.group2 <- data.frame(Gene = METABRIC.data[labels(MB1.metab.genes.cv.heat$rowDendrogram[[1]]),1], cv = metab.genes[match(c(METABRIC.data[labels(MB1.metab.genes.cv.heat$rowDendrogram[[1]]),1]),metab.genes$Gene),2], stringsAsFactors = FALSE)

# average expression of heatmap groups in individual samples

all.clinical.df2$MB1.bicluster.cv.heat.group1_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB1.bicluster.group1$Gene),-1]))
all.clinical.df2$MB1.bicluster.cv.heat.group2_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB1.bicluster.group2$Gene),-1]))
all.clinical.df2$MB1.all.cv.heat.group1_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB1.hi.cv.group1$Gene),-1]))
all.clinical.df2$MB1.all.cv.heat.group2_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB1.hi.cv.group2$Gene),-1]))
all.clinical.df2$MB1.metab.group1_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB1.metab.group1$Gene),-1]))
all.clinical.df2$MB1.metab.group2_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB1.metab.group2$Gene),-1]))


###MB3 (former Mito1): ----
# 1. get genes and their CVs in the bicluster:
MB3.bicluster.genes <- data.frame(Gene=(METABRIC.data[,1][mito.genes.loc]), cv=as.vector(mito.cv[[1]])[mito.genes.loc], stringsAsFactors = FALSE) 

# 2. get genes with the highest CV from the whole genome:
MB3.hi.cv.loc <- order(abs(mito.cv[[1]]),decreasing = TRUE)[seq(1000)]
MB3.hi.cv.genes <- data.frame(Gene = (METABRIC.data[,1][MB3.hi.cv.loc]), cv=as.vector(mito.cv[[1]])[MB3.hi.cv.loc], stringsAsFactors = FALSE) 

# 3. find the two anti-correlating groups by clustering the genes in the top 10 samples by heatmap.2
# a. for genes with the highest CV from the whole genome:

# pdf("MB3.all.cv.heat.pdf")

MB3.all.cv.heat <- heatmap.2(cor(t(METABRIC.data[,-1][MB3.hi.cv.loc,mito.sort[[1]][seq(10)]])),
                              trace = 'none')

# dev.off()

# note ICT1.sort[[1]][seq(10)] is ICT1.top.seed

# b. for the highest correlating genes in the bicluster (note that here instead of simple cor we use HclustGenesHiCor which trims the dendrogram to the highest correlating groups)



MB3.hicor.genes <- as.numeric(HclustGenesHiCor(METABRIC.data[,-1][mito.genes.loc,],
                                                mito.top.seed[[1]],
                                                cuts = 8))
# pdf("MB3.biclustergenes.cv.heat.pdf")

MB3.biclustergenes.cv.heat <- heatmap.2(cor(t(METABRIC.data[,-1][mito.genes.loc[MB3.hicor.genes],
                                                                  mito.sort[[1]][seq(10)]])),
                                         trace = 'none')

# dev.off()

# metabolic genes highly correlating with the bicluster 
MB3.hicor.metab.genes <- as.numeric(HclustGenesHiCor(METABRIC.data[,-1][metab.loc,],
                                                      mito.top.seed[[1]],
                                                      cuts = 6))

# pdf("MB3.metab.genes.cv.heat.pdf")

MB3.metab.genes.cv.heat <- heatmap.2(cor(t(METABRIC.data[,-1][metab.loc[MB3.hicor.metab.genes],
                                                               mito.sort[[1]][seq(10)]])),
                                      trace = 'none')
# dev.off()

# get genes together in dataframes with names and CVs in a list from all the above

df.list$MB3.hi.cv.genes.cv <- as.data.frame(MB3.hi.cv.genes) # genes with the highest CV from the whole genome:
df.list$MB3.bicluster.genes.cv <-  as.data.frame(MB3.bicluster.genes) # genes and their CVs in the bicluster
df.list$MB3.allgenome.cv <- data.frame(Gene = METABRIC.data[,1],CV = mito.cv[[1]]) #CVs of all genes in the genome
df.list$MB3.all.metab.genes.cv <- data.frame(Gene = as.character(metab.genes[,1]), CV = df.list$MB3.allgenome.cv[match(metab.genes[,1], df.list$MB3.allgenome.cv$Gene),2], stringsAsFactors = FALSE) #CVs of all metabolic genes in the genome
df.list$MB3.bicluster.group1 <- data.frame(Gene = METABRIC.data[labels(MB3.biclustergenes.cv.heat$rowDendrogram[[2]]),1], cv = MB3.bicluster.genes[match(c(METABRIC.data[labels(MB3.biclustergenes.cv.heat$rowDendrogram[[2]]),1]),MB3.bicluster.genes$Gene),2], stringsAsFactors = FALSE)
df.list$MB3.bicluster.group2 <- data.frame(Gene = METABRIC.data[labels(MB3.biclustergenes.cv.heat$rowDendrogram[[1]]),1], cv = MB3.bicluster.genes[match(c(METABRIC.data[labels(MB3.biclustergenes.cv.heat$rowDendrogram[[1]]),1]),MB3.bicluster.genes$Gene),2], stringsAsFactors = FALSE)
df.list$MB3.hi.cv.group1 <- data.frame(Gene = METABRIC.data[labels(MB3.all.cv.heat$rowDendrogram[[2]]),1], cv = MB3.hi.cv.genes[match(c(METABRIC.data[labels(MB3.all.cv.heat$rowDendrogram[[2]]),1]),MB3.hi.cv.genes$Gene),2], stringsAsFactors = FALSE)
df.list$MB3.hi.cv.group2 <- data.frame(Gene = METABRIC.data[labels(MB3.all.cv.heat$rowDendrogram[[1]]),1], cv = MB3.hi.cv.genes[match(c(METABRIC.data[labels(MB3.all.cv.heat$rowDendrogram[[1]]),1]),MB3.hi.cv.genes$Gene),2], stringsAsFactors = FALSE)
df.list$MB3.metab.group1 <- data.frame(Gene = METABRIC.data[labels(MB3.metab.genes.cv.heat$rowDendrogram[[2]]),1], cv = metab.genes[match(c(METABRIC.data[labels(MB3.metab.genes.cv.heat$rowDendrogram[[2]]),1]),metab.genes$Gene),2], stringsAsFactors = FALSE)
df.list$MB3.metab.group2 <- data.frame(Gene = METABRIC.data[labels(MB3.metab.genes.cv.heat$rowDendrogram[[1]]),1], cv = metab.genes[match(c(METABRIC.data[labels(MB3.metab.genes.cv.heat$rowDendrogram[[1]]),1]),metab.genes$Gene),2], stringsAsFactors = FALSE)

# average expression of heatmap groups in individual samples

all.clinical.df2$MB3.bicluster.cv.heat.group1_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB3.bicluster.group1$Gene),-1]))
all.clinical.df2$MB3.bicluster.cv.heat.group2_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB3.bicluster.group2$Gene),-1]))
all.clinical.df2$MB3.all.cv.heat.group1_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB3.hi.cv.group1$Gene),-1]))
all.clinical.df2$MB3.all.cv.heat.group2_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB3.hi.cv.group2$Gene),-1]))
all.clinical.df2$MB3.metab.group1_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB3.metab.group1$Gene),-1]))
all.clinical.df2$MB3.metab.group2_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB3.metab.group2$Gene),-1]))

###MB2 (former Mito2):----
# 1. get genes and their CVs in the bicluster:
MB2.bicluster.genes <- data.frame(Gene=(METABRIC.data[,1][mito.genes.loc]), cv=as.vector(mito.cv[[2]])[mito.genes.loc], stringsAsFactors = FALSE) 

# 2. get genes with the highest CV from the whole genome:
MB2.hi.cv.loc <- order(abs(mito.cv[[2]]),decreasing = TRUE)[seq(1000)]
MB2.hi.cv.genes <- data.frame(Gene = (METABRIC.data[,1][MB2.hi.cv.loc]), cv=as.vector(mito.cv[[2]])[MB2.hi.cv.loc], stringsAsFactors = FALSE) 

# 3. find the two anti-correlating groups by clustering the genes in the top 10 samples by heatmap.2
# a. for genes with the highest CV from the whole genome:

# pdf("MB2.all.cv.heat.pdf")

MB2.all.cv.heat <- heatmap.2(cor(t(METABRIC.data[,-1][MB2.hi.cv.loc,mito.sort[[2]][seq(10)]])),
                               trace = 'none')
# dev.off()

# note ICT1.sort[[1]][seq(10)] is ICT1.top.seed

# b. for the highest correlating genes in the bicluster (note that here instead of simple cor we use HclustGenesHiCor which trims the dendrogram to the highest correlating groups)

MB2.hicor.genes <- as.numeric(HclustGenesHiCor(METABRIC.data[,-1][mito.genes.loc,],
                                                 mito.top.seed[[2]],
                                                 cuts = 4))
# pdf("MB2.biclustergenes.cv.heat.pdf")

MB2.biclustergenes.cv.heat <- heatmap.2(cor(t(METABRIC.data[,-1][mito.genes.loc[MB2.hicor.genes],
                                                                   mito.sort[[2]][seq(10)]])),
                                          trace = 'none')

# dev.off()

# metabolic genes highly correlating with the bicluster 
MB2.hicor.metab.genes <- as.numeric(HclustGenesHiCor(METABRIC.data[,-1][metab.loc,],
                                                       mito.top.seed[[2]],
                                                       cuts = 6))

# pdf("MB2.metab.genes.cv.heat.pdf")

MB2.metab.genes.cv.heat <- heatmap.2(cor(t(METABRIC.data[,-1][metab.loc[MB2.hicor.metab.genes],
                                                                mito.sort[[2]][seq(10)]])),
                                       trace = 'none')

# dev.off()

# get genes together in dataframes with names and CVs in a list from all the above

df.list$MB2.hi.cv.genes.cv <- as.data.frame(MB2.hi.cv.genes) # genes with the highest CV from the whole genome:
df.list$MB2.bicluster.genes.cv <-  as.data.frame(MB2.bicluster.genes) # genes and their CVs in the bicluster
df.list$MB2.allgenome.cv <- data.frame(Gene = METABRIC.data[,1],CV = mito.cv[[2]]) #CVs of all genes in the genome
df.list$MB2.all.metab.genes.cv <- data.frame(Gene = as.character(metab.genes[,1]), CV = df.list$MB2.allgenome.cv[match(metab.genes[,1], df.list$MB2.allgenome.cv$Gene),2], stringsAsFactors = FALSE) #CVs of all metabolic genes in the genome
df.list$MB2.bicluster.group1 <- data.frame(Gene = METABRIC.data[labels(MB2.biclustergenes.cv.heat$rowDendrogram[[2]]),1], cv = MB2.bicluster.genes[match(c(METABRIC.data[labels(MB2.biclustergenes.cv.heat$rowDendrogram[[2]]),1]),MB2.bicluster.genes$Gene),2], stringsAsFactors = FALSE)
df.list$MB2.bicluster.group2 <- data.frame(Gene = METABRIC.data[labels(MB2.biclustergenes.cv.heat$rowDendrogram[[1]]),1], cv = MB2.bicluster.genes[match(c(METABRIC.data[labels(MB2.biclustergenes.cv.heat$rowDendrogram[[1]]),1]),MB2.bicluster.genes$Gene),2], stringsAsFactors = FALSE)
df.list$MB2.hi.cv.group1 <- data.frame(Gene = METABRIC.data[labels(MB2.all.cv.heat$rowDendrogram[[2]]),1], cv = MB2.hi.cv.genes[match(c(METABRIC.data[labels(MB2.all.cv.heat$rowDendrogram[[2]]),1]),MB2.hi.cv.genes$Gene),2], stringsAsFactors = FALSE)
df.list$MB2.hi.cv.group2 <- data.frame(Gene = METABRIC.data[labels(MB2.all.cv.heat$rowDendrogram[[1]]),1], cv = MB2.hi.cv.genes[match(c(METABRIC.data[labels(MB2.all.cv.heat$rowDendrogram[[1]]),1]),MB2.hi.cv.genes$Gene),2], stringsAsFactors = FALSE)
df.list$MB2.metab.group1 <- data.frame(Gene = METABRIC.data[labels(MB2.metab.genes.cv.heat$rowDendrogram[[2]]),1], cv = metab.genes[match(c(METABRIC.data[labels(MB2.metab.genes.cv.heat$rowDendrogram[[2]]),1]),metab.genes$Gene),2], stringsAsFactors = FALSE)
df.list$MB2.metab.group2 <- data.frame(Gene = METABRIC.data[labels(MB2.metab.genes.cv.heat$rowDendrogram[[1]]),1], cv = metab.genes[match(c(METABRIC.data[labels(MB2.metab.genes.cv.heat$rowDendrogram[[1]]),1]),metab.genes$Gene),2], stringsAsFactors = FALSE)




# average expression of heatmap groups in individual samples

all.clinical.df2$MB2.bicluster.cv.heat.group1_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB2.bicluster.group1$Gene),-1]))
all.clinical.df2$MB2.bicluster.cv.heat.group2_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB2.bicluster.group2$Gene),-1]))
all.clinical.df2$MB2.all.cv.heat.group1_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB2.hi.cv.group1$Gene),-1]))
all.clinical.df2$MB2.all.cv.heat.group2_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB2.hi.cv.group2$Gene),-1]))
all.clinical.df2$MB2.metab.group1_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB2.metab.group1$Gene),-1]))
all.clinical.df2$MB2.metab.group2_exp <- as.numeric(colMeans(METABRIC.data[which(METABRIC.data[,1] %in% df.list$MB2.metab.group2$Gene),-1]))


###### Prolifation Index ----
prol.genes <- c('BIRC5','CCNB1','CDC20','NUF2','CEP55','NDC80','MKI67',
                'PTTG1','RRM2','TYMS','UBE2C')
prol.loc <- which(METABRIC.data[,1] %in% prol.genes)
# note no PTTG1 measured in METABRIC data
all.clinical.df2$proliferation.index <- as.numeric(colMeans(METABRIC.data[prol.loc,-1]))

###### Claudin low ------


METABRIC_CL_Fougner.df <- read_excel('CL_CoreCL_METABRIC_Fougner.xlsx', sheet = 2 )
METABRIC_CL_Fougner.df2 <- METABRIC_CL_Fougner.df[,c(1,10,13,33:35)]

all.clinical.df2 <- left_join(all.clinical.df2, METABRIC_CL_Fougner.df2, by = c('sample' = 'Patient_ID'))

CL_1_3.df <-  read_excel("METABRIC_samples.xlsx")
all.clinical.df2 <- left_join(all.clinical.df2, CL_1_3.df, by = c('sample' = 'samples'))



####### Pommier developmental stages -----

Pommier_dev_genesets <- as.list(read_excel('Pommier_dev_genesets.xlsx'))


tmp1 <- as.matrix(METABRIC.data[,-1])
row.names(tmp1) <- METABRIC.data$external_gene_name

ssGSEA.test <- gsva(tmp1,
                    gset.idx.list = Pommier_dev_genesets,
                    method = 'gsva')

ssGSEA.test.df <- rownames_to_column(as.data.frame(t(ssGSEA.test)), var = 'samples' ) %>% mutate(samples = str_replace(samples, "MB.", "MB-"))

all.clinical.df2 <- left_join(all.clinical.df2, ssGSEA.test.df, by = c('sample' = 'samples'))

all.clinical.df2 <- all.clinical.df2 %>% mutate(MaSC_sum = (MaSC_up - MaSC_down)/2) %>% 
  mutate(LP_sum = (LP_up - LP_down)/2) %>% 
  mutate(mL_sum = (mL_up - mL_down)/2)

# save-----

df.list <- append(df.list,list("allSampleParameters"=all.clinical.df2), 0)

WriteXLS(df.list,ExcelFileName = 'METABRIC_allSampleParameters_allBiclustergenelists_corrected_groups2.xlsx')

rm(chr.loc, ssGSEA.test, i, attributes,Clinical_Overall_Survival_Data_from_METABRIC,Complete_METABRIC_Clinical_Features_Data, Complete_METABRIC_Clinical_Survival_Data_from_METABRIC,Complete_METABRIC_Copy_Number_Data,ensembl,METABRIC.copy,METABRIC.estimate, probe2gene,probe2gene2,sort.df,tmp.df,i,chr.loc,new.order,new.order.vec, unique.loc)

save.image(file="METABRIC_starting_data_final_corr_groups.Rdata")

save(df.list, file = 'essential_METABRIC_MCbiclust_data.Rdata')

