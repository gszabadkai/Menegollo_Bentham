# Code to look at PARADIGM results and their distribution among the bicluster forks

library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(WriteXLS)
library(pheatmap)
library(matrixStats)
library(ggstatsplot)

# 1. Load data ----
# PARADIGM datafile obtained from https://pubmed.ncbi.nlm.nih.gov/29622464/ - files available upon request from GS
# PanGyn_PARADIGM_2173.txt
paradigm.2173 <- read.table('PanGyn_PARADIGM_2173.txt',header = TRUE, sep = '\t')

# TCGA RNAseq data (all biclusters combined)
load('../TCGA.all.biclusters.RNAseq.Rdata')

# create df with only fork samples and most different paradigm values ----

fork.loc <- which(TCGA.all.biclusters.RNAseq.df2$MB1.fork != 'None' |
                    TCGA.all.biclusters.RNAseq.df2$MB3.fork != 'None' |
                    TCGA.all.biclusters.RNAseq.df2$MB2.fork != 'None')
fork.samps <- as.character(TCGA.all.biclusters.RNAseq.df2$X[fork.loc])
fork.samps2 <- sapply(strsplit(fork.samps,'-'), function(x) paste(x[seq(3)],collapse = '.'))
fork.pd.loc <- which(colnames(paradigm.2173) %in% fork.samps2)
paradigm.forks.2173 <- paradigm.2173[,fork.pd.loc]

# find paradigm IPLs with highest variability between all fork samples ----

paradigm.forks.2173.1K <-  order(rowVars(as.matrix(paradigm.forks.2173), na.rm = TRUE), decreasing = TRUE)[1:250]
paradigm.forks.2173.1K.df <- paradigm.forks.2173[paradigm.forks.2173.1K,]

save(paradigm.forks.2173.1K.df, file = "paradigm.forks.250.Rdata")

# create annotation from bicluster data for fork and all BRCA samples ----

TCGA.all.biclusters.RNAseq.df2.t <- t(TCGA.all.biclusters.RNAseq.df2)
samples <-  sapply(strsplit(TCGA.all.biclusters.RNAseq.df2.t['sample', ],'-'), function(x) paste(x[seq(3)],collapse = '.')) 
colnames(TCGA.all.biclusters.RNAseq.df2.t) <- samples
phenotype <- as.data.frame(t(TCGA.all.biclusters.RNAseq.df2.t[ ,colnames(paradigm.forks.2173.1K.df)]))

# convert dataframes to matrix for heatmap generation
paradigm.forks.2173.1K.matrix <-   as.matrix(paradigm.forks.2173.1K.df)

# create matrices for proliferation and tumour purity

prol.matrix <- t(as.matrix(as.numeric(phenotype[,35])))
colnames(prol.matrix) <- rownames(phenotype)

purity.matrix <- t(as.matrix(as.numeric(phenotype[,36])))
colnames(purity.matrix) <- rownames(phenotype)
purity.matrix[is.na(purity.matrix)] = 0.5

# Fig 3A create heatmap (only fork samples) ----

# annotation
fork_col <- structure(names = c("Upper", "Lower", 'None'), c('#bd0026','#253494', "ghostwhite"))
ha.forks = HeatmapAnnotation(
  MB1.fork = phenotype[[5]],
  MB2.fork = phenotype[[11]],
  MB3.fork = phenotype[[8]],
  PAM50 = phenotype[[12]],
  Histology = phenotype[[27]],
  col = list(MB1.fork = fork_col,
             MB2.fork = fork_col,
             MB3.fork = fork_col,
             PAM50 = structure(names = c("Basal", "LumA", "LumB", "Her2", "Normal"), 
                                      brewer.pal(5, "Set1")),
             Histology = structure(names = c("Invasive Ductal Carcinoma", "Invasive Lobular Carcinoma", "Mixed", "Other"), 
                                   brewer.pal(4, "Dark2"))),
  na_col = "grey", 
  border = TRUE,
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    MB1.fork = list(title = "MB1 forks"),
    MB2.fork = list(title = "MB2 forks"),
    MB3.fork = list(title = "MB3 forks"),
    PAM50 = list(title = "PAM50"),
    Histology = list(title = "Histology"))
)

# heatmap(s)
# distance and clustering methods were chosen to best reflect differences betwwen forks

col_fun = colorRamp2(c(-3, 0, 3), c("#377EB8", "white", "#E41A1C"))
ht_1_par = Heatmap(paradigm.forks.2173.1K.matrix,
                  # breaks = (seq(-5,5,length.out = 100)),
                  col = col_fun, 
                  name = "Paradigm",
                  clustering_distance_columns = "pearson",
                  clustering_distance_rows = "pearson",
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2",
                  show_row_dend = T, 
                  show_column_dend = T,
                  show_column_names = FALSE,
                  show_row_names = FALSE,
                  row_dend_side = "right",
                  show_heatmap_legend = F,
                  row_split = 11, row_title = "P%s", row_title_rot = 0, row_title_gp = gpar(col = "#377EB8",fontsize = 10, fontface= 'bold'),
                  # row_title = "Paradigm concepts with highest variance (n = 250, P1-11 clusters)",
                  column_split = 4,
                  # row_km = 12,
                  # column_km = 6,
                  row_gap = unit(1.5, "mm"),
                  # bottom_annotation = ha.forks, 
                  column_title = NULL
                  # row_split = factor(cgi2, levels = c("Island", "Shore", "Shelf", "OpenSea")), 
                  # row_title_gp = gpar(col = "#FFFFFF00"))
)

ht_2_prol <- Heatmap(prol.matrix, name = "Proliferation", 
                     col = colorRamp2(c(-1, 1), c("white", "#00441b")), height = unit(5, "mm"),
                     row_title = "Proliferation", 
                     row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                     row_title_rot = 0)

ht_3_purity <- Heatmap(purity.matrix, name = "Tumour purity", 
                       col = colorRamp2(c(0, 1), c("white", "#e6550d")), height = unit(5, "mm"),
                       row_title = "Tumour purity", 
                       row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                       row_title_rot = 0)
# draw(ht_3_purity)
ht_fork_list <- ht_1_par %v% ht_2_prol %v% ht_3_purity %v% ha.forks 
pdf("fork_ht2.pdf", width = 12, height = 10)
draw(ht_fork_list, annotation_legend_side = "bottom", heatmap_legend_side = "bottom", adjust_annotation_extension = FALSE)

annotation_titles = c(MB1.fork = 'MB1 forks',
                      MB2.fork = 'MB2 forks',
                      MB3.fork = 'MB3 forks',
                      PAM50 = 'PAM50',
                      Histology = "Histology")
for(an in names(annotation_titles)) {
  decorate_annotation(an, {
    grid.text(annotation_titles[an], unit(-2, "mm"), just = "right", gp = gpar(fontsize=10, fontface= 'bold'))
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}

decorate_annotation("PAM50", {
  grid.lines(unit(c(-30, 0), "mm"), unit(c(1, 1), "npc"))
})

decorate_annotation("MB1.fork", {
  grid.lines(unit(c(-30, 0), "mm"), unit(c(1, 1), "npc"))
})

dev.off()


# get paradigm names in different clusters

paradigm.clusters.fork.list <- list()
paradigm.clusters.fork.list$"slice1" <- rownames(paradigm.forks.2173.1K.df)[row_order(ht_1_par)[[1]]]
paradigm.clusters.fork.list$"slice2" <- rownames(paradigm.forks.2173.1K.df)[row_order(ht_1_par)[[2]]]
paradigm.clusters.fork.list$"slice3" <- rownames(paradigm.forks.2173.1K.df)[row_order(ht_1_par)[[3]]]
paradigm.clusters.fork.list$"slice4" <- rownames(paradigm.forks.2173.1K.df)[row_order(ht_1_par)[[4]]]
paradigm.clusters.fork.list$"slice5" <- rownames(paradigm.forks.2173.1K.df)[row_order(ht_1_par)[[5]]]
paradigm.clusters.fork.list$"slice6" <- rownames(paradigm.forks.2173.1K.df)[row_order(ht_1_par)[[6]]]
paradigm.clusters.fork.list$"slice7" <- rownames(paradigm.forks.2173.1K.df)[row_order(ht_1_par)[[7]]]
paradigm.clusters.fork.list$"slice8" <- rownames(paradigm.forks.2173.1K.df)[row_order(ht_1_par)[[8]]]
paradigm.clusters.fork.list$"slice9" <- rownames(paradigm.forks.2173.1K.df)[row_order(ht_1_par)[[9]]]
paradigm.clusters.fork.list$"slice10" <- rownames(paradigm.forks.2173.1K.df)[row_order(ht_1_par)[[10]]]
paradigm.clusters.fork.list$"slice11" <- rownames(paradigm.forks.2173.1K.df)[row_order(ht_1_par)[[11]]]
  

# Fig S3 ggstatplot

# summarise parameters in clusters, compare their means between forks

paradigm.forks.2173.1K.df2 <- paradigm.forks.2173.1K.df %>% rownames_to_column('P') 

paradigm.forks.2173.1K.df2  <- mutate(paradigm.forks.2173.1K.df2, Pcluster = ifelse(P %in% paradigm.clusters.fork.list$slice1, 'P1',
                                                                                    ifelse(P %in% paradigm.clusters.fork.list$slice2, 'P2',
                                                                                           ifelse(P %in% paradigm.clusters.fork.list$slice3, 'P3',
                                                                                                  ifelse(P %in% paradigm.clusters.fork.list$slice4, 'P4',
                                                                                                         ifelse(P %in% paradigm.clusters.fork.list$slice5, 'P5',
                                                                                                                ifelse(P %in% paradigm.clusters.fork.list$slice6, 'P6',
                                                                                                                       ifelse(P %in% paradigm.clusters.fork.list$slice7, 'P7',
                                                                                                                              ifelse(P %in% paradigm.clusters.fork.list$slice8, 'P8',
                                                                                                                                     ifelse(P %in% paradigm.clusters.fork.list$slice9, 'P9',
                                                                                                                                            ifelse(P %in% paradigm.clusters.fork.list$slice10, 'P10',
                                                                                                                                                   ifelse(P %in% paradigm.clusters.fork.list$slice11, 'P11', 'none'
))))))))))))

temp <- paradigm.forks.2173.1K.df2[,c(1,414)]

WriteXLS(temp, ExcelFileName = 'paradigm_clusters.xls')



paradigm.forks.2173.1K.df2.sum_by_P <- as.data.frame(t(paradigm.forks.2173.1K.df2 %>% group_by(Pcluster) %>% 
  summarise_at(vars(c(2:413)), funs(mean)) %>% 
  column_to_rownames(var = 'Pcluster')))
paradigm.forks.2173.1K.df2.sum_by_P <- paradigm.forks.2173.1K.df2.sum_by_P[,c(1,4:11,2,3)]


# compare means of clustermeans/sample between groups
         

paradigm.forks.2173.1K.df2.sum_by_P <- rownames_to_column(paradigm.forks.2173.1K.df2.sum_by_P)
phenotype <- rownames_to_column(phenotype)


paradigm.statplot.df <- left_join(phenotype, paradigm.forks.2173.1K.df2.sum_by_P)

paradigm.statplot.df$MB1.forkscale.fork_0.4 <- paste("MB1", paradigm.statplot.df$MB1.forkscale.fork_0.4,sep="_")
paradigm.statplot.df$MB2.forkscale.fork_0.4 <- paste("MB2", paradigm.statplot.df$MB2.forkscale.fork_0.4,sep="_")
paradigm.statplot.df$MB3.forkscale.fork_0.4 <- paste("MB3", paradigm.statplot.df$MB3.forkscale.fork_0.4,sep="_")

paradigm.statplot.df[,c(202:207)] <- lapply(paradigm.statplot.df[,c(202:207)], function(x) {gsub("Lower", "LF", x)})
paradigm.statplot.df[,c(202:207)] <- lapply(paradigm.statplot.df[,c(202:207)], function(x) {gsub("Upper", "UF", x)})

save(paradigm.statplot.df, file='TCGA.all.biclusters.phenotypes.w.paradigm.cluster.values.Rdata')

paradigm.statplot.df.g <- gather(paradigm.statplot.df, 'MB1.forkscale.fork_0.4', 'MB2.forkscale.fork_0.4', 'MB3.forkscale.fork_0.4', key = 'forkscales', value = 'forks')

paradigm.statplot.df.g2 <- gather(paradigm.statplot.df.g, P1:P11, key = 'Pcluster', value = 'Pvalue')

pdf('..', width = 40, height = 27)

grouped_ggbetweenstats(paradigm.statplot.df.g2 %>% subset(forks %in% c('MB1_UF','MB1_LF', 'MB2_LF','MB2_UF', 'MB3_LF','MB3_UF')), 
               x = forks, 
               y = Pvalue,
               grouping.var = Pcluster,
               type = 'np')

dev.off()
