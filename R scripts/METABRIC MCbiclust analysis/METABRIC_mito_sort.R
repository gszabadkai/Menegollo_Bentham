args<-commandArgs(TRUE)

library(MCbiclust)

load('../METABRIC_DATA.RData')
load('..METABRIC_mito_p3.RData')

pattern <- as.numeric(args[1])
pat.tot <- 3

av.corvec <- average.corvec[[pattern]]
top.genes.num <- 750
top.genes <- order(abs(av.corvec),decreasing = T)[seq(length = top.genes.num)]
top.mat <- as.matrix(METABRIC.data[,-1])[top.genes,]

top.sort <- SampleSort(top.mat, seed = top.seed[[pattern]])

save(list = c('top.sort'),
     file = paste('METABRIC_mito_sort',
                  '_pattern',pattern,
                  '_of',pat.tot,'.RData',sep = ''))