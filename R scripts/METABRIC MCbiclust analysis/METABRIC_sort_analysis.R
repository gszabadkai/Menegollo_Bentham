library(MCbiclust)

# Load data
load('../METABRIC_DATA.RData')
load('../METABRIC_gene_lists.RData')
load('../METABRIC_sort_data.RData')

# PC1 calculation
mito.pc1 <- list()
for(i in 1:3){
  top.mat <- METABRIC.data[order(abs(mito.cv[[i]]),decreasing = TRUE)[seq(1000)],-1]
  mito.pc1[[i]] <- PC1VecFun(top.mat,seed.sort = mito.sort[[i]],n = 10)}


random.pc1 <- list()
for(i in 1:2){
  top.mat <- METABRIC.data[order(abs(random.cv[[i]]),decreasing = TRUE)[seq(1000)],-1]
  random.pc1[[i]] <- PC1VecFun(top.mat,seed.sort = random.sort[[i]],n = 10)}

ICT1.pc1 <- list()
for(i in 1:1){
  top.mat <- METABRIC.data[order(abs(ICT1.cv[[i]]),decreasing = TRUE)[seq(1000)],-1]
  ICT1.pc1[[i]] <- PC1VecFun(top.mat,seed.sort = ICT1.sort[[i]],n = 10)}


# CV comparison

all.cv <- data.frame(mito.cv, ICT1.cv, random.cv)
colnames(all.cv) <- c('Mito1','Mito2','Mito3','ICT1','Rand1','Rand2')

library(gdata)
mito.genes <- as.character(read.xls('Data/Human.MitoCarta2.0.xls', sheet = 2)[,4])
mito.genes.loc <- which(METABRIC.data[,1] %in% mito.genes)

CVPlot(cv.df = all.cv, geneset.loc = mito.genes.loc,
       geneset.name = 'Mito',
       cnames =  c('Mito1','Mito2','Mito3','ICT1','Rand1','Rand2'))

# ICT1 ~ Mito2
# Rand1 = Mito1
# Rand2 = Mito3

# Gene set enrichment

all.gsea <- list()
for(i in 1:6){
all.gsea[[i]] <- GOEnrichmentAnalysis(gene.names = as.character(METABRIC.data[,1]),
                     gene.values = all.cv[,i],sig.rate = 0.05)
}


save(list = c('mito.pc1','random.pc1','ICT1.pc1','all.gsea'),
     file = 'METABRIC_PC1_GSEA.RData')


write.table(all.gsea[[1]],col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = '\t',file = 'RData/mito1_gsea.txt')
write.table(all.gsea[[2]],col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = '\t',file = 'RData/mito2_gsea.txt')
write.table(all.gsea[[3]],col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = '\t',file = 'RData/mito3_gsea.txt')
write.table(all.gsea[[4]],col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = '\t',file = 'RData/ICT1_gsea.txt')
