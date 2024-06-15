library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)

# Load data
load("../METABRIC_starting_data_final_corr_groups.Rdata")

# Select 7K top and bottom genes by MB1 CV

MB1.hi.cv.loc.14K <- c(order(ICT1.cv[[1]], decreasing = TRUE)[seq(7000)], rev(order(ICT1.cv[[1]], decreasing = FALSE)[seq(7000)]))

# set up gene lists 
mitocarta2.0 <- read_excel('Human.MitoCarta2.0.xls', sheet=2)
mito.loc <- which(METABRIC.data[,1] %in% mitocarta2.0[,4])

metab.genes.xls <- read_excel('GS metabolic genes list.xlsx', sheet=1)
metab.loc <- which(METABRIC.data[,1] %in% metab.genes.xls[,2])

# Create Complexheatmap Fig 1G

# for annotation
phenotype2.CVs <-  data_frame('genes' = METABRIC.data[,1],
                             'MB1.cv' = ICT1.cv[[1]],
                             'MB2.cv' = mito.cv[[2]],
                             'MB3.cv' = mito.cv[[1]])

# divide positive and negative CV genes in the bicluster genes

MB1.pos.bicluster.genes <- MB1.bicluster.genes[MB1.bicluster.genes$cv > 0,]
MB1.neg.bicluster.genes <- MB1.bicluster.genes[MB1.bicluster.genes$cv < 0,]
MB2.pos.bicluster.genes <- MB2.bicluster.genes[MB2.bicluster.genes$cv > 0,]
MB2.neg.bicluster.genes <- MB2.bicluster.genes[MB2.bicluster.genes$cv < 0,]
MB3.pos.bicluster.genes <- MB3.bicluster.genes[MB3.bicluster.genes$cv > 0,]
MB3.neg.bicluster.genes <- MB3.bicluster.genes[MB3.bicluster.genes$cv < 0,]

phenotype2.CVs <-  phenotype2.CVs %>%
   mutate('Mitochondrial_genes' = ifelse(phenotype2.CVs$genes %in% mitocarta2.0[,4], 'Y', 'N')) %>%
   mutate('Metabolic_genes' = ifelse(phenotype2.CVs$genes %in% metab.genes.xls[,2], 'Y', 'N')) %>%
   mutate('MB1_bicluster_genes' = ifelse(phenotype2.CVs$genes %in% MB1.pos.bicluster.genes$Gene, 'U', ifelse(phenotype2.CVs$genes %in% MB1.neg.bicluster.genes$Gene, 'L', 'N'))) %>% 
   mutate('MB2_bicluster_genes' = ifelse(phenotype2.CVs$genes %in% MB2.pos.bicluster.genes$Gene, 'U', ifelse(phenotype2.CVs$genes %in% MB2.neg.bicluster.genes$Gene, 'L', 'N'))) %>% 
   mutate('MB3_bicluster_genes' = ifelse(phenotype2.CVs$genes %in% MB3.pos.bicluster.genes$Gene, 'U', ifelse(phenotype2.CVs$genes %in% MB3.neg.bicluster.genes$Gene, 'L', 'N')))

phenotype2.CVs <- phenotype2.CVs[MB1.hi.cv.loc.14K,]

#  color scheme
cat_col_metab_2 <- structure(names = c("Y", "N"), c("#238b45", "#f0f0f0"))
cat_col_mito_2 <- structure(names = c("Y", "N"), c("#ff7f00", "#f0f0f0"))
cat_col_MB1_2 <- structure(names = c("U", "L", "N"), c('#D55E00', '#a6bddb',"white"))
cat_col_MB2_2 <- structure(names = c("U", "L", "N"), c('#1b7837', '#980043',"white"))
cat_col_MB3_2 <- structure(names = c("U", "L", "N"), c('#762a83', '#d7301f',"white"))

#  annotations:
ha.bottom <-  HeatmapAnnotation(
   MB1_biclust = phenotype2.CVs$MB1_bicluster_genes,
   MB1.CV = anno_points(phenotype2.CVs$MB1.cv,
                        border = FALSE,
                        gp = gpar(col = '#737373', alpha = 0.1, size = 0.1),
                        height = unit(1.5, "cm"),
                        axis_param = list(side = "right", at = c(-1, 0, 1))),
   MB2_biclust = phenotype2.CVs$MB2_bicluster_genes,
   MB2.CV = anno_points(phenotype2.CVs$MB2.cv,
                        border = FALSE,
                        gp = gpar(col = '#737373', alpha = 0.1, size = 0.1),
                        height = unit(1.5, "cm"),
                        axis_param = list(side = "right", at = c(-1, 0, 1))),
   MB3_biclust = phenotype2.CVs$MB3_bicluster_genes,
   MB3.CV = anno_points(phenotype2.CVs$MB3.cv,
                        border = FALSE,
                        gp = gpar(col = '#737373', alpha = 0.1, size = 0.1),
                        height = unit(1.5, "cm"),
                        axis_param = list(side = "right", at = c(-1, 0, 1))),
   col = list(MB1_biclust = cat_col_MB1_2, MB2_biclust = cat_col_MB2_2, MB3_biclust = cat_col_MB3_2),
   show_annotation_name = FALSE,
   gap = unit(2, "mm")
)

ha.top <- HeatmapAnnotation(
   Metabolic = phenotype2.CVs$Metabolic_genes,
   Mitochondrial = phenotype2.CVs$Mitochondrial_genes,
   col = list(Mitochondrial = cat_col_mito_2, Metabolic = cat_col_metab_2),
   show_legend = c(F,F),
   border = FALSE,
   show_annotation_name = FALSE,
   gap = unit(1.5, "mm"))

col_fun = colorRamp2(c(-0.8,-0.3, 0.3, 0.8), c("#253494","white","white","#bd0026"))

# heatmap
ht <- Heatmap(cor(t(METABRIC.data[,-1][MB1.hi.cv.loc.14K,ICT1.sort[[1]][seq(10)]])),
              name = 'Pearson r',
              cluster_rows = F,
              cluster_columns = F,
              col = col_fun,
              row_split = c(rep('LF', 7000),rep('UF', 7000)), row_gap = unit(2, "mm"),
              column_split = c(rep('LF', 7000),rep('UF', 7000)), column_gap = unit(2, "mm"),
              # row_title = c('MB1_UF','MB1_LF'),
              # column_title = c('MB1_UF','MB1_LF'),
              show_row_names = FALSE,
              show_column_names = FALSE,
              bottom_annotation = ha.bottom,
              top_annotation = ha.top,
)
pdf("..", height = 8, width = 12)

draw(ht, padding = unit(c(2, 2, 2, 40), "mm"), annotation_legend_side = "right")

dev.off()

# the heatmap is further edited in Adobe illustrator

# ______________________________________________________________________

# Fig. 1H
# Create data to assess mitochondrial gene enrichment in forks 
# create lists of genes for overall intersections: intersections between hi.cv, biclusters, metab genes and mitocarta2 in each bicluster----

venn.sets <- list()

venn.sets$wg = as.character(df.list$MB1.allgenome.cv[,1])
venn.sets$wg.mito.genes = as.character(mitocarta2.0[,4])
venn.sets$wg.metab.genes = metab.genes[,1]
venn.sets$wg.mito.metab.genes = intersect(as.character(mitocarta2.0[,4]),metab.genes[,1])
venn.sets$wg.non.mito.metab.genes = setdiff(union(metab.genes[,1],mitocarta2.0[,4]), mitocarta2.0[,4])
venn.sets$wg.mito.non.metab.genes = setdiff(union(metab.genes[,1],mitocarta2.0[,4]), metab.genes[,1])

venn.sets$MB1.group1 = df.list$MB1.hi.cv.group1[,1]
venn.sets$MB1.group1.mito.genes = intersect(df.list$MB1.hi.cv.group1[,1],mitocarta2.0[,4])
venn.sets$MB1.group1.metab.genes = intersect(df.list$MB1.hi.cv.group1[,1],metab.genes[,1])
venn.sets$MB1.group1.mito.metab.genes = intersect(df.list$MB1.hi.cv.group1[,1], intersect(as.character(mitocarta2.0[,4]),metab.genes[,1]))
venn.sets$MB1.group1.non.mito.metab.genes = setdiff(venn.sets$MB1.group1.metab.genes, venn.sets$MB1.group1.mito.metab.genes)
venn.sets$MB1.group1.mito.non.metab.genes = setdiff(venn.sets$MB1.group1.mito.genes, venn.sets$MB1.group1.mito.metab.genes)

venn.sets$MB1.group2 = df.list$MB1.hi.cv.group2[,1]
venn.sets$MB1.group2.mito.genes = intersect(df.list$MB1.hi.cv.group2[,1],mitocarta2.0[,4])
venn.sets$MB1.group2.metab.genes = intersect(df.list$MB1.hi.cv.group2[,1],metab.genes[,1])
venn.sets$MB1.group2.mito.metab.genes = intersect(df.list$MB1.hi.cv.group2[,1], intersect(as.character(mitocarta2.0[,4]),metab.genes[,1]))
venn.sets$MB1.group2.non.mito.metab.genes = setdiff(venn.sets$MB1.group2.metab.genes, venn.sets$MB1.group2.mito.metab.genes)
venn.sets$MB1.group2.mito.non.metab.genes = setdiff(venn.sets$MB1.group2.mito.genes, venn.sets$MB1.group2.mito.metab.genes)

venn.sets$MB2.group1 = df.list$MB2.hi.cv.group1[,1]
venn.sets$MB2.group1.mito.genes = intersect(df.list$MB2.hi.cv.group1[,1],mitocarta2.0[,4])
venn.sets$MB2.group1.metab.genes = intersect(df.list$MB2.hi.cv.group1[,1],metab.genes[,1])
venn.sets$MB2.group1.mito.metab.genes = intersect(df.list$MB2.hi.cv.group1[,1], intersect(as.character(mitocarta2.0[,4]),metab.genes[,1]))
venn.sets$MB2.group1.non.mito.metab.genes = setdiff(venn.sets$MB2.group1.metab.genes, venn.sets$MB2.group1.mito.metab.genes)
venn.sets$MB2.group1.mito.non.metab.genes = setdiff(venn.sets$MB2.group1.mito.genes, venn.sets$MB2.group1.mito.metab.genes)

venn.sets$MB2.group2 = df.list$MB2.hi.cv.group2[,1]
venn.sets$MB2.group2.mito.genes = intersect(df.list$MB2.hi.cv.group2[,1],mitocarta2.0[,4])
venn.sets$MB2.group2.metab.genes = intersect(df.list$MB2.hi.cv.group2[,1],metab.genes[,1])
venn.sets$MB2.group2.mito.metab.genes = intersect(df.list$MB2.hi.cv.group2[,1], intersect(as.character(mitocarta2.0[,4]),metab.genes[,1]))
venn.sets$MB2.group2.non.mito.metab.genes = setdiff(venn.sets$MB2.group2.metab.genes, venn.sets$MB2.group2.mito.metab.genes)
venn.sets$MB2.group2.mito.non.metab.genes = setdiff(venn.sets$MB2.group2.mito.genes, venn.sets$MB2.group2.mito.metab.genes)

venn.sets$MB3.group1 = df.list$MB3.hi.cv.group1[,1]
venn.sets$MB3.group1.mito.genes = intersect(df.list$MB3.hi.cv.group1[,1],mitocarta2.0[,4])
venn.sets$MB3.group1.metab.genes = intersect(df.list$MB3.hi.cv.group1[,1],metab.genes[,1])
venn.sets$MB3.group1.mito.metab.genes = intersect(df.list$MB3.hi.cv.group1[,1], intersect(as.character(mitocarta2.0[,4]),metab.genes[,1]))
venn.sets$MB3.group1.non.mito.metab.genes = setdiff(venn.sets$MB3.group1.metab.genes, venn.sets$MB3.group1.mito.metab.genes)
venn.sets$MB3.group1.mito.non.metab.genes = setdiff(venn.sets$MB3.group1.mito.genes, venn.sets$MB3.group1.mito.metab.genes)

venn.sets$MB3.group2 = df.list$MB3.hi.cv.group2[,1]
venn.sets$MB3.group2.mito.genes = intersect(df.list$MB3.hi.cv.group2[,1],mitocarta2.0[,4])
venn.sets$MB3.group2.metab.genes = intersect(df.list$MB3.hi.cv.group2[,1],metab.genes[,1])
venn.sets$MB3.group2.mito.metab.genes = intersect(df.list$MB3.hi.cv.group2[,1], intersect(as.character(mitocarta2.0[,4]),metab.genes[,1]))
venn.sets$MB3.group2.non.mito.metab.genes = setdiff(venn.sets$MB3.group2.metab.genes, venn.sets$MB3.group2.mito.metab.genes)
venn.sets$MB3.group2.mito.non.metab.genes = setdiff(venn.sets$MB3.group2.mito.genes, venn.sets$MB3.group2.mito.metab.genes)

# calculate the relative enrichment of mitochondrial and metabolic genes in the high CV gene sets, show it in heatmap ----
# create matrix of gene numbers in each group
sets <- t(matrix(sapply(venn.sets, FUN=length),6,7))

# calculate 'enrichment'
norm.sets <- sets[,-1] / sets[,1]
norm.sets.2 <- t(t(norm.sets)[,-1] / t(norm.sets)[,1])

# reorder
norm.sets.2[ , c(2,5)] <- norm.sets.2[ , c(5,2)]
norm.sets.2[ , c(1,2)] <- norm.sets.2[ , c(2,1)]
norm.sets.2[ , c(4,5)] <- norm.sets.2[ , c(5,4)]

# create table with numbers for the heatmap
sets.2 <- sets[-1,-1]
sets.2[ , c(2,5)] <- sets.2[ , c(5,2)]
sets.2[ , c(1,2)] <- sets.2[ , c(2,1)]
sets.2[ , c(4,5)] <- sets.2[ , c(5,4)]

# create heatmap

pdf("..", height = 8, width = 16)

norm.sets.2.heatmap <- heatmap.2(norm.sets.2,
          Rowv = F,
          Colv = F,
          dendrogram = "none",
          trace="none",
          rowsep = c(2,4),
          colsep = c(1,2,3,4),
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          margins = c(2,12),
          col = colorRampPalette(c("#003366", "ghostwhite", "#990000"))(n = 299),
          breaks = c(seq(0, 0.7,length=100),seq(0.71, 1.3,length=100),seq(1.31,3.5,length=100)),
          cellnote = sets.2,
          notecol = "black",
          notecex = 1.9,
          keysize = 1.2,
          key.title = "",
          key.xlab = "",
          key.ylab = "",
          density.info="none",
          lhei = c(1.1,2),
          labRow = "",
          labCol = "")
dev.off()


# the heatmap was further annotated in illustrator


