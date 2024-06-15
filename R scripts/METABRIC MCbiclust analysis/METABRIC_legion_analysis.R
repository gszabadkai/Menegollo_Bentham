load('../METABRIC_DATA.RData')
load('../METABRIC_gene_lists.RData')

X <- as.matrix(METABRIC.data[,-1])

f1 <- 'METABRIC_random_'
f2 <- list.files()[grep(f1,list.files())]

all.cor.vec <- list()
all.seed <- list()
all.hicor <- list()
for(i in seq_len(length(f2))){
  load(f2[i])
  all.cor.vec[[i]] <- test.cor.vec
  all.seed[[i]] <- test.seed
  all.hicor[[i]] <- test.hicor.genes
}

all.cor.vec.mat <- as.matrix(as.data.frame(all.cor.vec))
rm(test.cor.vec, test.seed, test.hicor.genes, all.cor.vec)

clust.groups <- SilhouetteClustGroups(cor.vec.mat = all.cor.vec.mat,
                                      max.clusters = 20,
                                      plots = T, rand.vec = T)

for(i in seq_len(length(clust.groups))){
  group.max <- max(cor(as.matrix(all.cor.vec.mat[,clust.groups[[i]]])) -
                     diag(length(which(clust.groups[[i]] == TRUE))))
  
  group.min <- min(cor(as.matrix(all.cor.vec.mat[,clust.groups[[i]]])))
  
  minmax.diff <- group.max - group.min
  if(minmax.diff > 1.5){
    kmeans.test <- kmeans(t(all.cor.vec.mat[,clust.groups[[i]]]),
                          centers = 2)$cluster
    all.cor.vec.mat[,clust.groups[[i]]][,which(kmeans.test == 1)] <-
      - all.cor.vec.mat[,clust.groups[[i]]][,which(kmeans.test == 1)]
  }
}

av.cv.fun <- function(x){
  if(length(which(x== TRUE)) > 1){
    return(rowMeans(all.cor.vec.mat[,x]))
  }else{
    return(all.cor.vec.mat[,x])
  }
}

average.corvec <- lapply(clust.groups, FUN = av.cv.fun)

top.genes.num <- 750
top.seed <- list()
for(pattern in seq_len(length(clust.groups))){
  av.corvec <- average.corvec[[pattern]]
  groups <- clust.groups[[pattern]]
  
  top.genes <- order(abs(av.corvec),decreasing = T)[seq(length = top.genes.num)]
  top.seed.score <- seq(length = length(which(groups == TRUE)))
  
  for(i in seq(length = length(top.seed.score))){
    l1 <- top.genes
    l2 <- all.seed[groups][[i]]
    top.seed.score[i] <- mean(abs(cor(t(as.matrix(X)[l1,l2]))))
  }
  top.seed[[pattern]] <- all.seed[groups][[which.max(top.seed.score)[1]]]
}

save(list = c('top.seed','average.corvec'),
     file = paste(f1,'p',length(clust.groups),
                  '.RData',sep = ''))





