args<-commandArgs(TRUE)

library(MCbiclust)

load('../METABRIC_DATA.RData')
load('../METABRIC_gene_lists.RData')

run.seed <- as.numeric(args[1])

set.seed(run.seed)

X <- as.matrix(METABRIC.data[mito.genes.loc,-1])

test.seed <- FindSeed(gem = X, seed.size = 10,
                      iterations = 2000)

test.hicor.genes <- as.numeric(HclustGenesHiCor(gem = X,
                                                seed = test.seed,
                                                cuts = 8))

test.top.mat <- X[test.hicor.genes,]

test.cor.vec <- CVEval(gem.part = test.top.mat,
                       gem.all = as.matrix(METABRIC.data[,-1]),
                       seed = test.seed,
                       splits = 8)

save(list = c('test.cor.vec','test.seed','test.hicor.genes'),
     file = paste('METABRIC_mito_seed_',run.seed,'.RData',sep = ''))

