# Calculate pathway activities to compare pathways between forks and Pam50 categories

library(ComplexHeatmap)
library(tidyr)
library(dplyr)
library(tidyverse)
library(stringr)
library(RColorBrewer)

# Load data
load("../METABRIC_starting_data_final_corr_groups.Rdata")

# glutamine metabolism ----
# create df with expression of a pathway in each sample and the PAM50 and fork categories
glut.df <- (subset(METABRIC.data, external_gene_name %in% c('GFPT1', 'GFPT2', 'PPAT', 'GLUL', 'GLS', 'GLS2', 'GLUD1', 'GLUD2', 'GAD1', 'GAD2', 'GCLC', 'GCLM', 'GSS', 'GGCT', 'GPT', 'GPT2', 'PSAT1')))
glut.df2 <- as.data.frame(t(glut.df[,-1]))
colnames(glut.df2) <- glut.df[,1]
glut.df2 <- rownames_to_column(glut.df2, var = 'sample')
glut.df2$PAM50 <- all.clinical.df2$PAM50
glut.df2$MB1 <- all.clinical.df2$MB1.forkscale.fork_0.8


# combine categories and calculate mean by each category for each enzyme expression
glut.df2.means <- glut.df2 %>%unite(combined, PAM50, MB1, remove = FALSE) %>% 
  pivot_longer(c(PAM50, MB1, combined), names_to = "cat_key", values_to = "cat") %>% 
  group_by(cat) %>% 
  summarise_at(vars(GAD1:PSAT1), mean, na.rm = TRUE) %>% 
  column_to_rownames(var = 'cat')

glut.df2.sums <- glut.df2 %>%unite(combined, PAM50, MB1, remove = FALSE) %>% 
  pivot_longer(c(PAM50, MB1, combined), names_to = "cat_key", values_to = "cat") %>% 
  group_by(cat) %>% 
  summarise_at(vars(GAD1:PSAT1), sum, na.rm = TRUE) %>% 
  column_to_rownames(var = 'cat')

# order genes by CVs

cv.metab <- data.frame(df.list$MB1.all.metab.genes.cv, df.list$MB2.all.metab.genes.cv[,2], df.list$MB3.all.metab.genes.cv[,2], stringsAsFactors = FALSE)
colnames(cv.metab) <- c('Gene', 'MB1_CV', 'MB2_CV', 'MB3_CV')

glut.cv <- subset(cv.metab, Gene %in% c('GFPT1', 'GFPT2', 'PPAT', 'GLUL', 'GLS', 'GLS2', 'GLUD1', 'GLUD2', 'GAD1', 'GAD2', 'GCLC', 'GCLM', 'GSS', 'GGCT', 'GPT', 'GPT2', 'PSAT1'))
glut.cv.MB1 <- glut.cv[order(-glut.cv$MB1_CV),]

glut.df2.means <- glut.df2.means[glut.cv.MB1$Gene]
glut.df2.sums <- glut.df2.sums[glut.cv.MB1$Gene]


#  for heatmap: compare MB1_UF, MB1_LF with Pam50 categories

glut.df2.means.heatmap <- glut.df2.means[c('MB1_UF', 'MB1_LF', "LumA", "LumA_MB1_LF", "LumA_MB1_none", "LumB", "LumB_MB1_UF", "LumB_MB1_none"),]
glut.df2.sums.heatmap <- glut.df2.sums[c('MB1_UF', 'MB1_LF', "LumA", "LumA_MB1_LF", "LumA_MB1_none", "LumB", "LumB_MB1_UF", "LumB_MB1_none"),]

# ComplexHeatmap with CVs

row_ha = rowAnnotation("MB1\nCVs" = anno_barplot(glut.cv.MB1$MB1_CV, border = T, bar_width = 0.5, annotation_name_gp = gpar(fontface = "bold")))

Heatmap(t(glut.df2.means.heatmap),
        name = "METABRIC \nnormalised \ngene expression",
        row_order = colnames(glut.df2.means.heatmap),
        row_names_side = "left",
        width = unit(4, "cm"), height = unit(8, "cm"),
        column_km = 3, 
        column_title = c("LF", "UF", "PAM50"),
        right_annotation = row_ha
        
        # column_names_gp = gpar(col = c("green", "orange", "purple"))
        )

# including all groups
glut.df2.allgroups <- glut.df2

glut.df2.allgroups$MB2 <- all.clinical.df2$MB2.forkscale.fork_0.8
glut.df2.allgroups$MB3 <- all.clinical.df2$MB3.forkscale.fork_0.8

glut.df2.means.allgroups <- glut.df2.allgroups %>% unite(combined, PAM50, MB1, MB2, MB3, sep = " - ", remove = FALSE) %>% 
  pivot_longer(c(PAM50, MB1, MB2, MB3, combined), names_to = "cat_key", values_to = "cat") %>% 
  group_by(cat) %>% 
  summarise_at(vars(GAD1:PSAT1), mean, na.rm = TRUE) %>% 
  column_to_rownames(var = 'cat')

Heatmap(glut.df2.means.allgroups,
        name = "METABRIC \nnormalised \ngene expression",
        width = unit(10, "cm"), height = unit(20, "cm")
)

# make the dataframe only with all combinations and without the NC group

glut.df2.means.allgroups.allcombinations <- glut.df2.means.allgroups %>% 
  rownames_to_column(var = "cat") %>% 
  filter(nchar(as.character(cat)) > 10) %>% 
  filter(!str_detect(cat, "NC"))

Heatmap(glut.df2.means.allgroups.allcombinations %>% column_to_rownames(var = 'cat'),
        name = "METABRIC \nnormalised \ngene expression",
        width = unit(10, "cm"), height = unit(20, "cm"),
        row_order = glut.df2.means.allgroups.allcombinations$cat
)

# create separate columns for PAM50 and the three biclusters for annotation of the heatmap


glut.df2.means.allgroups.allcombinations.annotations <- glut.df2.means.allgroups.allcombinations %>% 
  mutate('PAM50' = ifelse(grepl('Basal', `cat`), 'Basal',
                          ifelse(grepl('LumA', `cat`), 'LumA',
                                 ifelse(grepl('LumB', `cat`), 'LumB',
                                        ifelse(grepl('Her2', `cat`), 'Her2',
                                               ifelse(grepl('Normal', `cat`), 'Normal', 'NC')))))) %>% 
  mutate('MB1' = ifelse(grepl('MB1_UF', `cat`), 'UF',
                         ifelse(grepl('MB1_LF', `cat`), 'LF',
                                ifelse(grepl('MB1_none', `cat`), 'none','NC')))) %>% 
  mutate('MB2' = ifelse(grepl('MB2_UF', `cat`), 'UF',
                        ifelse(grepl('MB2_LF', `cat`), 'LF',
                               ifelse(grepl('MB2_none', `cat`), 'none','NC')))) %>% 
  mutate('MB3' = ifelse(grepl('MB3_UF', `cat`), 'UF',
                        ifelse(grepl('MB3_LF', `cat`), 'LF',
                               ifelse(grepl('MB3_none', `cat`), 'none','NC')))) %>% 
  column_to_rownames(var = 'cat')

glut.annotations <- glut.df2.means.allgroups.allcombinations.annotations[,c("PAM50","MB1","MB2","MB3")]

# colors used in all heatmaps ----


fork_col <- structure(names = c("UF", "LF", 'none'), c('#bd0026','#253494', "ghostwhite"))
PAM50.col <- structure(names = c("Basal", "LumA", "LumB", "Her2", "Normal"), brewer.pal(5, "Set1"))

# glutamine heatmaps ----

main <- as.matrix(glut.df2.means.allgroups.allcombinations %>% column_to_rownames(var = 'cat'))

# the following heatmap can be plotted either with row clustering (this revelas grouping by biclusters) or by the default rwo order
# by Pam50 categories (this shows that Pam50 gives no information on metabolic patterns)
set.seed(123)

ht_list = Heatmap(main,
                  name = "METABRIC \nnormalised \ngene expression",
                  width = unit(10, "cm"), height = unit(12, "cm"),
                  # row_order = oxphos.df2.means.allgroups.allcombinations$cat,
                  column_names_gp = gpar(fontsize = 6),
                  clustering_distance_columns = "pearson",
                  clustering_distance_rows = "pearson",
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2",
                  row_names_side = "left",
                  column_names_side = "top",
                  show_row_names = FALSE,
                  show_column_names = TRUE,
                  row_title = NULL,
                  column_title = NULL,
                  column_split = 3,
                  # column_order = colnames(oxphos.df2.means.heatmap)
                  row_km = 3
                  ) +
  Heatmap(glut.df2.means.allgroups.allcombinations.annotations["PAM50"],
          name = "PAM50",
          width = unit(0.3, "cm"), height = unit(20, "cm"),
          col = PAM50.col,
          column_names_gp = gpar(fontsize = 9),
          show_row_names = FALSE,
          column_names_side = "top") +
  Heatmap(glut.df2.means.allgroups.allcombinations.annotations[,c("MB1","MB2","MB3")],
          name = "biclusters\nforks",
          width = unit(0.9, "cm"), height = unit(20, "cm"),
          col = fork_col,
          column_names_gp = gpar(fontsize = 9),
          show_row_names = FALSE,
          column_names_side = "top",
          cluster_columns = FALSE)

pdf("...", width = 8, height = 10)

draw(ht_list)  

dev.off()

# OXPHOS----

TEMP <- scan(text= 'COX4I2 NDUFC2 NDUFC1 NDUFB10 NDUFB9 NDUFB8 NDUFB7 NDUFB6 NDUFB5 NDUFB4 NDUFB3 NDUFB2 NDUFB1 NDUFA11 NDUFAB1 NDUFA10 NDUFA9 NDUFA8 NDUFA7 NDUFA6 NDUFA5 NDUFA4L2 NDUFA3 NDUFA2 NDUFA1 NDUFV3 NDUFV2 NDUFV1 NDUFS8 NDUFS7 NDUFS6 NDUFS5 NDUFS4 NDUFS3 NDUFS2 NDUFS1 ND6 ND5 ND4L ND4 ND3 ND2 ND1 COX17 COX15 COX11 COX8A COX7C COX7B COX7A1 COX6C CYC1 COX5B ATP5A1 ATP12A PPA1 SDHC SDHD SDHA SDHB UQCR11 UQCR10 UQCRQ UQCRB UQCRH UQCRC2 UQCRC1 CYTB UQCRFS1 COX6B1 COX6A1 COX5A COX2 COX10 COX1 COX3 ATP5B ATP5C1 ATP5O ATP5D ATP5E ATP5G2 ATP6 ATP5F1 ATP5I ATP5J ATP5J2 ATP8 ATP5H ATP5L ATP6V1A ATP6V1B1 ATP6V1C1 ATP6V1D ATP6V1E1 ATP6V1F ATP6V1G2 ATP6V1H ATP6V0A1 ATP6V0D1 ATP6V0E1 ATP6AP1 ATP6V0C NDUFB11 NDUFA12 NDUFA13 CYCS' , what="")


# create df with expression of a pathway in each sample and the PAM50 and fork categories
oxphos.df <- (subset(METABRIC.data, external_gene_name %in% TEMP))
oxphos.df2 <- as.data.frame(t(oxphos.df[,-1]))
colnames(oxphos.df2) <- oxphos.df[,1]
oxphos.df2 <- rownames_to_column(oxphos.df2, var = 'sample')
oxphos.df2$PAM50 <- all.clinical.df2$PAM50
oxphos.df2$MB1 <- all.clinical.df2$MB1.forkscale.fork_0.8


# combine categories and calculate mean by each category for each enzyme expression
oxphos.df2.means <- oxphos.df2 %>%unite(combined, PAM50, MB1, remove = FALSE) %>% 
  pivot_longer(c(PAM50, MB1, combined), names_to = "cat_key", values_to = "cat") %>% 
  group_by(cat) %>% 
  summarise_at(vars(ATP5A1:UQCRQ), mean, na.rm = TRUE) %>% 
  column_to_rownames(var = 'cat')

# order genes by CVs

cv.metab <- data.frame(df.list$MB1.all.metab.genes.cv, df.list$MB2.all.metab.genes.cv[,2], df.list$MB3.all.metab.genes.cv[,2], stringsAsFactors = FALSE)
colnames(cv.metab) <- c('Gene', 'MB1_CV', 'MB2_CV', 'MB3_CV')

oxphos.cv <- subset(cv.metab, Gene %in% TEMP)
oxphos.cv.MB1 <- oxphos.cv[order(-oxphos.cv$MB1_CV),]

oxphos.df2.means <- oxphos.df2.means[oxphos.cv.MB1$Gene]
# oxphos.df2.sums <- oxphos.df2.sums[oxphos.cv.MB1$Gene]


#  for heatmap: compare MB1_UF, MB1_LF with Pam50 categories

oxphos.df2.means.heatmap <- oxphos.df2.means[c('MB1_UF', 'MB1_LF', "LumA", "LumA_MB1_LF", "LumA_MB1_none", "LumB", "LumB_MB1_UF", "LumB_MB1_none"),]
# oxphos.df2.sums.heatmap <- oxphos.df2.sums[c('MB1_UF', 'MB1_LF', "LumA", "LumA_MB1_LF", "LumA_MB1_none", "LumB", "LumB_MB1_UF", "LumB_MB1_none"),]


# heatmap(t(oxphos.df2.means.heatmap), Rowv = NA, scale = c("column"))
# heatmap(t(oxphos.df2.sums.heatmap), Rowv = NA, scale = c("column"))

# ComplexHeatmap with CVs

row_ha = rowAnnotation("MB1\nCVs" = anno_barplot(oxphos.cv.MB1$MB1_CV, border = T, bar_width = 0.5, annotation_name_gp = gpar(fontface = "bold")))

Heatmap(t(oxphos.df2.means.heatmap),
        name = "METABRIC \nnormalised \ngene expression",
        row_order = colnames(oxphos.df2.means.heatmap),
        row_names_side = "left",
        width = unit(4, "cm"), height = unit(8, "cm"),
        column_km = 3, 
        column_title = c("LF", "PAM50", "UF"),
        right_annotation = row_ha,
        show_row_names = FALSE
        
        # column_names_gp = gpar(col = c("green", "orange", "purple"))
)

# including all groups
oxphos.df2.allgroups <- oxphos.df2

oxphos.df2.allgroups$MB2 <- all.clinical.df2$MB2.forkscale.fork_0.8
oxphos.df2.allgroups$MB3 <- all.clinical.df2$MB3.forkscale.fork_0.8

oxphos.df2.means.allgroups <- oxphos.df2.allgroups %>%unite(combined, PAM50, MB1, MB2, MB3, sep = " - ", remove = FALSE) %>% 
  pivot_longer(c(PAM50, MB1, MB2, MB3, combined), names_to = "cat_key", values_to = "cat") %>% 
  group_by(cat) %>% 
  summarise_at(vars(ATP5A1:UQCRQ), mean, na.rm = TRUE) %>% 
  column_to_rownames(var = 'cat')

Heatmap(oxphos.df2.means.allgroups,
        name = "METABRIC \nnormalised \ngene expression",
        width = unit(10, "cm"), height = unit(20, "cm")
)

# make the dataframe only with all combinations and without the NC group

oxphos.df2.means.allgroups.allcombinations <- oxphos.df2.means.allgroups %>% 
  rownames_to_column(var = "cat") %>% 
  filter(nchar(as.character(cat)) > 10) %>% 
  filter(!str_detect(cat, "NC"))

Heatmap(oxphos.df2.means.allgroups.allcombinations %>% column_to_rownames(var = 'cat'),
        name = "METABRIC \nnormalised \ngene expression",
        width = unit(10, "cm"), height = unit(20, "cm"),
        row_order = oxphos.df2.means.allgroups.allcombinations$cat
)

# create separate columns for PAM50 and the three biclusters for annotation of the heatmap


oxphos.df2.means.allgroups.allcombinations.annotations <- oxphos.df2.means.allgroups.allcombinations %>% 
  mutate('PAM50' = ifelse(grepl('Basal', `cat`), 'Basal',
                          ifelse(grepl('LumA', `cat`), 'LumA',
                                 ifelse(grepl('LumB', `cat`), 'LumB',
                                        ifelse(grepl('Her2', `cat`), 'Her2',
                                               ifelse(grepl('Normal', `cat`), 'Normal', 'NC')))))) %>% 
  mutate('MB1' = ifelse(grepl('MB1_UF', `cat`), 'UF',
                        ifelse(grepl('MB1_LF', `cat`), 'LF',
                               ifelse(grepl('MB1_none', `cat`), 'none','NC')))) %>% 
  mutate('MB2' = ifelse(grepl('MB2_UF', `cat`), 'UF',
                        ifelse(grepl('MB2_LF', `cat`), 'LF',
                               ifelse(grepl('MB2_none', `cat`), 'none','NC')))) %>% 
  mutate('MB3' = ifelse(grepl('MB3_UF', `cat`), 'UF',
                        ifelse(grepl('MB3_LF', `cat`), 'LF',
                               ifelse(grepl('MB3_none', `cat`), 'none','NC')))) %>% 
  column_to_rownames(var = 'cat')

# oxphos.annotations <- oxphos.df2.means.allgroups.allcombinations.annotations[,c("PAM50","MB1","MB2","MB3")]


fork_col <- structure(names = c("UF", "LF", 'none'), c('#bd0026','#253494', "ghostwhite"))
PAM50.col <- structure(names = c("Basal", "LumA", "LumB", "Her2", "Normal"), brewer.pal(5, "Set1"))


oxphos.main <- as.matrix(oxphos.df2.means.allgroups.allcombinations %>% column_to_rownames(var = 'cat'))

# oxphos heatmaps ----
# the following heatmap can be plotted either with row clustering (this revelas grouping by biclusters) or by the default rwo order
# by Pam50 categories (this shows that Pam50 gives no information on metabolic patterns)

ht_list.oxphos = Heatmap(oxphos.main,
                  name = "METABRIC \nnormalised \ngene expression",
                  width = unit(10, "cm"), height = unit(12, "cm"),
                  # row_order = oxphos.df2.means.allgroups.allcombinations$cat,
                  clustering_distance_columns = "pearson",
                  clustering_distance_rows = "pearson",
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2",
                  row_names_side = "left",
                  column_names_side = "top",
                  show_row_names = FALSE,
                  show_column_names = FALSE,
                  row_title = NULL,
                  column_title = NULL,
                  column_split = 4,
                  # column_order = colnames(oxphos.df2.means.heatmap)
                  row_km = 3
                                ) +
  Heatmap(oxphos.df2.means.allgroups.allcombinations.annotations["PAM50"],
          name = "PAM50",
          width = unit(0.3, "cm"), height = unit(20, "cm"),
          col = PAM50.col,
          column_names_gp = gpar(fontsize = 9),
          show_row_names = FALSE,
          column_names_side = "top") +
  Heatmap(oxphos.df2.means.allgroups.allcombinations.annotations[,c("MB1","MB2","MB3")],
          name = "biclusters\nforks",
          width = unit(0.9, "cm"), height = unit(20, "cm"),
          col = fork_col,
          show_row_names = FALSE,
          column_names_gp = gpar(fontsize = 9),
          column_names_side = "top",
          cluster_columns = FALSE)

pdf("..", width = 8, height = 10)
ht_oxphos <- draw(ht_list.oxphos)
ht_oxphos
dev.off()

column_order(ht_oxphos)


oxphos.clusters.fork.list <- list()
oxphos.clusters.fork.list$"slice1" <- colnames(oxphos.main)[column_order(ht_oxphos)[[1]]]
oxphos.clusters.fork.list$"slice2" <- colnames(oxphos.main)[column_order(ht_oxphos)[[2]]]
oxphos.clusters.fork.list$"slice3" <- colnames(oxphos.main)[column_order(ht_oxphos)[[3]]]
oxphos.clusters.fork.list$"slice4" <- colnames(oxphos.main)[column_order(ht_oxphos)[[4]]]


# pyruvate metabolism----

# create df with expression of a pathway in each sample and the PAM50 and fork categories
pyr.df <- (subset(METABRIC.data, external_gene_name %in% c('PCX','PCK2','MDH1','ME2','ME3','MOD1','ME1','ME3','PDHA1','PDHA2','PDHB','PDHX','DLAT','MDH1','MDH2','PC','PKLR','PKM2','PCK1','PCK2','LDHA','LDHAL6B','LDHB','LDHC','LDHD','ACSS1','ACSS2','ACACA','ACACB','ACAT1','ACOT12','ALDH1B1','ALDH2','ALDH3A2','ALDH7A1','ALDH9A1')))
pyr.df2 <- as.data.frame(t(pyr.df[,-1]))
colnames(pyr.df2) <- pyr.df[,1]
pyr.df2 <- rownames_to_column(pyr.df2, var = 'sample')
pyr.df2$PAM50 <- all.clinical.df2$PAM50
pyr.df2$MB1 <- all.clinical.df2$MB1.forkscale.fork_0.8


# combine categories and calculate mean by each category for each enzyme expression
pyr.df2.means <- pyr.df2 %>%unite(combined, PAM50, MB1, remove = FALSE) %>% 
  pivot_longer(c(PAM50, MB1, combined), names_to = "cat_key", values_to = "cat") %>% 
  group_by(cat) %>% 
  summarise_at(vars(ACACA:PKLR), mean, na.rm = TRUE) %>% 
  column_to_rownames(var = 'cat')

# pyr.df2.sums <- pyr.df2 %>%unite(combined, PAM50, MB1, remove = FALSE) %>% 
#   pivot_longer(c(PAM50, MB1, combined), names_to = "cat_key", values_to = "cat") %>% 
#   group_by(cat) %>% 
#   summarise_at(vars(GAD1:PSAT1), sum, na.rm = TRUE) %>% 
#   column_to_rownames(var = 'cat')

# order genes by CVs

cv.metab <- data.frame(df.list$MB1.all.metab.genes.cv, df.list$MB2.all.metab.genes.cv[,2], df.list$MB3.all.metab.genes.cv[,2], stringsAsFactors = FALSE)
colnames(cv.metab) <- c('Gene', 'MB1_CV', 'MB2_CV', 'MB3_CV')

pyr.cv <- subset(cv.metab, Gene %in% c('PCX','PCK2','MDH1','ME2','ME3','MOD1','ME1','ME3','PDHA1','PDHA2','PDHB','PDHX','DLAT','MDH1','MDH2','PC','PKLR','PKM2','PCK1','PCK2','LDHA','LDHAL6B','LDHB','LDHC','LDHD','ACSS1','ACSS2','ACACA','ACACB','ACAT1','ACOT12','ALDH1B1','ALDH2','ALDH3A2','ALDH7A1','ALDH9A1'))
pyr.cv.MB1 <- pyr.cv[order(-pyr.cv$MB1_CV),]

pyr.df2.means <- pyr.df2.means[pyr.cv.MB1$Gene]
# pyr.df2.sums <- pyr.df2.sums[pyr.cv.MB1$Gene]


#  for heatmap: compare MB1_UF, MB1_LF with Pam50 categories

pyr.df2.means.heatmap <- pyr.df2.means[c('MB1_UF', 'MB1_LF', "LumA", "LumA_MB1_LF", "LumA_MB1_none", "LumB", "LumB_MB1_UF", "LumB_MB1_none"),]
# pyr.df2.sums.heatmap <- pyr.df2.sums[c('MB1_UF', 'MB1_LF', "LumA", "LumA_MB1_LF", "LumA_MB1_none", "LumB", "LumB_MB1_UF", "LumB_MB1_none"),]


# heatmap(t(pyr.df2.means.heatmap), Rowv = NA, scale = c("column"))
# heatmap(t(pyr.df2.sums.heatmap), Rowv = NA, scale = c("column"))

# ComplexHeatmap with CVs

row_ha = rowAnnotation("MB1\nCVs" = anno_barplot(pyr.cv.MB1$MB1_CV, border = T, bar_width = 0.5, annotation_name_gp = gpar(fontface = "bold")))

Heatmap(t(pyr.df2.means.heatmap),
        name = "METABRIC \nnormalised \ngene expression",
        row_order = colnames(pyr.df2.means.heatmap),
        row_names_side = "left",
        width = unit(4, "cm"), height = unit(12, "cm"),
        column_km = 3, 
        column_title = c("UF", "PAM50", "LF"),
        right_annotation = row_ha,
        show_row_names = TRUE
        
        # column_names_gp = gpar(col = c("green", "orange", "purple"))
)

# including all groups
pyr.df2.allgroups <- pyr.df2

pyr.df2.allgroups$MB2 <- all.clinical.df2$MB2.forkscale.fork_0.8
pyr.df2.allgroups$MB3 <- all.clinical.df2$MB3.forkscale.fork_0.8

pyr.df2.means.allgroups <- pyr.df2.allgroups %>%unite(combined, PAM50, MB1, MB2, MB3, sep = " - ", remove = FALSE) %>% 
  pivot_longer(c(PAM50, MB1, MB2, MB3, combined), names_to = "cat_key", values_to = "cat") %>% 
  group_by(cat) %>% 
  summarise_at(vars(ACACA:PKLR), mean, na.rm = TRUE) %>% 
  column_to_rownames(var = 'cat')

# make the dataframe only with all combinations and without the NC group

pyr.df2.means.allgroups.allcombinations <- pyr.df2.means.allgroups %>% 
  rownames_to_column(var = "cat") %>% 
  filter(nchar(as.character(cat)) > 10) %>% 
  filter(!str_detect(cat, "NC"))


# create separate columns for PAM50 and the three biclusters for annotation of the heatmap


pyr.df2.means.allgroups.allcombinations.annotations <- pyr.df2.means.allgroups.allcombinations %>% 
  mutate('PAM50' = ifelse(grepl('Basal', `cat`), 'Basal',
                          ifelse(grepl('LumA', `cat`), 'LumA',
                                 ifelse(grepl('LumB', `cat`), 'LumB',
                                        ifelse(grepl('Her2', `cat`), 'Her2',
                                               ifelse(grepl('Normal', `cat`), 'Normal', 'NC')))))) %>% 
  mutate('MB1' = ifelse(grepl('MB1_UF', `cat`), 'UF',
                        ifelse(grepl('MB1_LF', `cat`), 'LF',
                               ifelse(grepl('MB1_none', `cat`), 'none','NC')))) %>% 
  mutate('MB2' = ifelse(grepl('MB2_UF', `cat`), 'UF',
                        ifelse(grepl('MB2_LF', `cat`), 'LF',
                               ifelse(grepl('MB2_none', `cat`), 'none','NC')))) %>% 
  mutate('MB3' = ifelse(grepl('MB3_UF', `cat`), 'UF',
                        ifelse(grepl('MB3_LF', `cat`), 'LF',
                               ifelse(grepl('MB3_none', `cat`), 'none','NC')))) %>% 
  column_to_rownames(var = 'cat')

# pyr.annotations <- pyr.df2.means.allgroups.allcombinations.annotations[,c("PAM50","MB1","MB2","MB3")]


fork_col <- structure(names = c("UF", "LF", 'none'), c('#bd0026','#253494', "ghostwhite"))
PAM50.col <- structure(names = c("Basal", "LumA", "LumB", "Her2", "Normal"), brewer.pal(5, "Set1"))


pyr.main <- as.matrix(pyr.df2.means.allgroups.allcombinations %>% column_to_rownames(var = 'cat'))

# pyr heatmaps ----
# the following heatmap can be plotted either with row clustering (this revelas grouping by biclusters) or by the default rwo order
# by Pam50 categories (this shows that Pam50 gives no information on metabolic patterns)

ht_list.pyr = Heatmap(pyr.main,
                         name = "METABRIC \nnormalised \ngene expression",
                         width = unit(10, "cm"), height = unit(20, "cm"),
                         row_order = pyr.df2.means.allgroups.allcombinations$cat,
                         clustering_distance_columns = "pearson",
                         clustering_distance_rows = "pearson",
                         # clustering_method_rows = "ward.D2",
                         # clustering_method_columns = "ward.D2",
                         row_names_side = "left",
                         column_names_side = "top",
                         show_row_names = FALSE,
                         show_column_names = TRUE,
                         row_title = NULL,
                         column_title = NULL,
                         column_split = 3,
                         # column_order = colnames(pyr.df2.means.heatmap)
                         # row_km = 4
) +
  Heatmap(pyr.df2.means.allgroups.allcombinations.annotations["PAM50"],
          name = "PAM50",
          width = unit(10/17, "cm"), height = unit(20, "cm"),
          col = PAM50.col,
          show_row_names = FALSE,
          column_names_side = "top") +
  Heatmap(pyr.df2.means.allgroups.allcombinations.annotations[,c("MB1","MB2","MB3")],
          name = "biclusters\nforks",
          width = unit(30/17, "cm"), height = unit(20, "cm"),
          col = fork_col,
          show_row_names = FALSE,
          column_names_side = "top",
          cluster_columns = FALSE)


draw(ht_list.pyr)  

ht_pyr <- Heatmap(pyr.main,
                     name = "METABRIC \nnormalised \ngene expression",
                     width = unit(10, "cm"), height = unit(20, "cm"),
                     # row_order = pyr.df2.means.allgroups.allcombinations$cat,
                     clustering_distance_columns = "pearson",
                     clustering_distance_rows = "pearson",
                     clustering_method_rows = "ward.D2",
                     clustering_method_columns = "ward.D2",
                     row_names_side = "left",
                     column_names_side = "top",
                     show_row_names = FALSE,
                     show_column_names = FALSE,
                     row_title = NULL,
                     column_title = NULL,
                     column_split = 4,
                     # column_order = colnames(pyr.df2.means.heatmap)
                     row_km = 3)

column_order(ht_pyr)


pyr.clusters.fork.list <- list()
pyr.clusters.fork.list$"slice1" <- colnames(pyr.main)[column_order(ht_pyr)[[1]]]
pyr.clusters.fork.list$"slice2" <- colnames(pyr.main)[column_order(ht_pyr)[[2]]]
pyr.clusters.fork.list$"slice3" <- colnames(pyr.main)[column_order(ht_pyr)[[3]]]
pyr.clusters.fork.list$"slice4" <- colnames(pyr.main)[column_order(ht_pyr)[[4]]]





# TCA----

# create df with expression of a pathway in each sample and the PAM50 and fork categories
TCA.df <- (subset(METABRIC.data, external_gene_name %in% c('PDHA1','PDHA2','PDHB','PDHX','CS','ACO1','ACO2','IDH1','IDH2','IDH3A','IDH3B','IDH3G','OGDH','OGDHL','DLST','SUCLG1','SUCLG2','SDHA','SDHB','SDHC','SDHD','FH1','DHTKD1','L2HGDH','DLD','MDH1','MDH2','PCK1','PCK2','ACLY','DLAT','PC')))
TCA.df2 <- as.data.frame(t(TCA.df[,-1]))
colnames(TCA.df2) <- TCA.df[,1]
TCA.df2 <- rownames_to_column(TCA.df2, var = 'sample')
TCA.df2$PAM50 <- all.clinical.df2$PAM50
TCA.df2$MB1 <- all.clinical.df2$MB1.forkscale.fork_0.8


# combine categories and calculate mean by each category for each enzyme expression
TCA.df2.means <- TCA.df2 %>%unite(combined, PAM50, MB1, remove = FALSE) %>% 
  pivot_longer(c(PAM50, MB1, combined), names_to = "cat_key", values_to = "cat") %>% 
  group_by(cat) %>% 
  summarise_at(vars(ACLY:SUCLG2), mean, na.rm = TRUE) %>% 
  column_to_rownames(var = 'cat')

# order genes by CVs

cv.metab <- data.frame(df.list$MB1.all.metab.genes.cv, df.list$MB2.all.metab.genes.cv[,2], df.list$MB3.all.metab.genes.cv[,2], stringsAsFactors = FALSE)
colnames(cv.metab) <- c('Gene', 'MB1_CV', 'MB2_CV', 'MB3_CV')

TCA.cv <- subset(cv.metab, Gene %in% c('PDHA1','PDHA2','PDHB','PDHX','CS','ACO1','ACO2','IDH1','IDH2','IDH3A','IDH3B','IDH3G','OGDH','OGDHL','DLST','SUCLG1','SUCLG2','SDHA','SDHB','SDHC','SDHD','FH1','DHTKD1','L2HGDH','DLD','MDH1','MDH2','PCK1','PCK2','ACLY','DLAT','PC'))
TCA.cv.MB1 <- TCA.cv[order(-TCA.cv$MB1_CV),]

TCA.df2.means <- TCA.df2.means[TCA.cv.MB1$Gene]
# TCA.df2.sums <- TCA.df2.sums[TCA.cv.MB1$Gene]


#  for heatmap: compare MB1_UF, MB1_LF with Pam50 categories

TCA.df2.means.heatmap <- TCA.df2.means[c('MB1_UF', 'MB1_LF', "LumA", "LumA_MB1_LF", "LumA_MB1_none", "LumB", "LumB_MB1_UF", "LumB_MB1_none"),]
# TCA.df2.sums.heatmap <- TCA.df2.sums[c('MB1_UF', 'MB1_LF', "LumA", "LumA_MB1_LF", "LumA_MB1_none", "LumB", "LumB_MB1_UF", "LumB_MB1_none"),]


# heatmap(t(TCA.df2.means.heatmap), Rowv = NA, scale = c("column"))
# heatmap(t(TCA.df2.sums.heatmap), Rowv = NA, scale = c("column"))

# ComplexHeatmap with CVs

row_ha = rowAnnotation("MB1\nCVs" = anno_barplot(TCA.cv.MB1$MB1_CV, border = T, bar_width = 0.5, annotation_name_gp = gpar(fontface = "bold")))

Heatmap(t(TCA.df2.means.heatmap),
        name = "METABRIC \nnormalised \ngene expression",
        row_order = colnames(TCA.df2.means.heatmap),
        row_names_side = "left",
        width = unit(4, "cm"), height = unit(12, "cm"),
        column_km = 3, 
        column_title = c("LF", "UF", "PAM50"),
        right_annotation = row_ha,
        show_row_names = TRUE
        
        # column_names_gp = gpar(col = c("green", "orange", "purple"))
)

# including all groups
TCA.df2.allgroups <- TCA.df2

TCA.df2.allgroups$MB2 <- all.clinical.df2$MB2.forkscale.fork_0.8
TCA.df2.allgroups$MB3 <- all.clinical.df2$MB3.forkscale.fork_0.8

TCA.df2.means.allgroups <- TCA.df2.allgroups %>%unite(combined, PAM50, MB1, MB2, MB3, sep = " - ", remove = FALSE) %>% 
  pivot_longer(c(PAM50, MB1, MB2, MB3, combined), names_to = "cat_key", values_to = "cat") %>% 
  group_by(cat) %>% 
  summarise_at(vars(ACLY:SUCLG2), mean, na.rm = TRUE) %>% 
  column_to_rownames(var = 'cat')

# make the dataframe only with all combinations and without the NC group

TCA.df2.means.allgroups.allcombinations <- TCA.df2.means.allgroups %>% 
  rownames_to_column(var = "cat") %>% 
  filter(nchar(as.character(cat)) > 10) %>% 
  filter(!str_detect(cat, "NC"))

# create separate columns for PAM50 and the three biclusters for annotation of the heatmap


TCA.df2.means.allgroups.allcombinations.annotations <- TCA.df2.means.allgroups.allcombinations %>% 
  mutate('PAM50' = ifelse(grepl('Basal', `cat`), 'Basal',
                          ifelse(grepl('LumA', `cat`), 'LumA',
                                 ifelse(grepl('LumB', `cat`), 'LumB',
                                        ifelse(grepl('Her2', `cat`), 'Her2',
                                               ifelse(grepl('Normal', `cat`), 'Normal', 'NC')))))) %>% 
  mutate('MB1' = ifelse(grepl('MB1_UF', `cat`), 'UF',
                        ifelse(grepl('MB1_LF', `cat`), 'LF',
                               ifelse(grepl('MB1_none', `cat`), 'none','NC')))) %>% 
  mutate('MB2' = ifelse(grepl('MB2_UF', `cat`), 'UF',
                        ifelse(grepl('MB2_LF', `cat`), 'LF',
                               ifelse(grepl('MB2_none', `cat`), 'none','NC')))) %>% 
  mutate('MB3' = ifelse(grepl('MB3_UF', `cat`), 'UF',
                        ifelse(grepl('MB3_LF', `cat`), 'LF',
                               ifelse(grepl('MB3_none', `cat`), 'none','NC')))) %>% 
  column_to_rownames(var = 'cat')

# TCA.annotations <- TCA.df2.means.allgroups.allcombinations.annotations[,c("PAM50","MB1","MB2","MB3")]


fork_col <- structure(names = c("UF", "LF", 'none'), c('#bd0026','#253494', "ghostwhite"))
PAM50.col <- structure(names = c("Basal", "LumA", "LumB", "Her2", "Normal"), brewer.pal(5, "Set1"))


TCA.main <- as.matrix(TCA.df2.means.allgroups.allcombinations %>% column_to_rownames(var = 'cat'))

# TCA heatmaps ----
# the following heatmap can be plotted either with row clustering (this revelas grouping by biclusters) or by the default rwo order
# by Pam50 categories (this shows that Pam50 gives no information on metabolic patterns)

ht_list.TCA = Heatmap(TCA.main,
                      name = "METABRIC \nnormalised \ngene expression",
                      width = unit(10, "cm"), height = unit(12, "cm"),
                      # row_order = TCA.df2.means.allgroups.allcombinations$cat,
                      column_names_gp = gpar(fontsize = 6),
                      clustering_distance_columns = "pearson",
                      clustering_distance_rows = "pearson",
                      clustering_method_rows = "ward.D2",
                      clustering_method_columns = "ward.D2",
                      row_names_side = "left",
                      column_names_side = "top",
                      show_row_names = FALSE,
                      show_column_names = TRUE,
                      row_title = NULL,
                      column_title = NULL,
                      column_split = 3,
                      # column_order = colnames(TCA.df2.means.heatmap)
                      row_km = 3
) +
  Heatmap(TCA.df2.means.allgroups.allcombinations.annotations["PAM50"],
          name = "PAM50",
          width = unit(0.3, "cm"), height = unit(20, "cm"),
          col = PAM50.col,
          show_row_names = FALSE,
          column_names_gp = gpar(fontsize = 9),
          column_names_side = "top") +
  Heatmap(TCA.df2.means.allgroups.allcombinations.annotations[,c("MB1","MB2","MB3")],
          name = "biclusters\nforks",
          width = unit(0.9, "cm"), height = unit(20, "cm"),
          col = fork_col,
          show_row_names = FALSE,
          column_names_gp = gpar(fontsize = 9),
          column_names_side = "top",
          cluster_columns = FALSE)

pdf("..", width = 8, height = 10)

draw(ht_list.TCA) 

dev.off()

ht_TCA <- Heatmap(TCA.main,
                  name = "METABRIC \nnormalised \ngene expression",
                  width = unit(10, "cm"), height = unit(20, "cm"),
                  # row_order = TCA.df2.means.allgroups.allcombinations$cat,
                  clustering_distance_columns = "pearson",
                  clustering_distance_rows = "pearson",
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2",
                  row_names_side = "left",
                  column_names_side = "top",
                  show_row_names = FALSE,
                  show_column_names = FALSE,
                  row_title = NULL,
                  column_title = NULL,
                  column_split = 4,
                  # column_order = colnames(TCA.df2.means.heatmap)
                  row_km = 3)

column_order(ht_TCA)


TCA.clusters.fork.list <- list()
TCA.clusters.fork.list$"slice1" <- colnames(TCA.main)[column_order(ht_TCA)[[1]]]
TCA.clusters.fork.list$"slice2" <- colnames(TCA.main)[column_order(ht_TCA)[[2]]]
TCA.clusters.fork.list$"slice3" <- colnames(TCA.main)[column_order(ht_TCA)[[3]]]
TCA.clusters.fork.list$"slice4" <- colnames(TCA.main)[column_order(ht_TCA)[[4]]]



# glycolysis----

# create df with expression of a pathway in each sample and the PAM50 and fork categories
glyc.df <- (subset(METABRIC.data, external_gene_name %in% c('GCK', 'HK1', 'HK2', 'HK3', 'GPI1', 'PFKL', 'PFKM', 'PFKP', 'PFKFB1', 'PFKFB2', 'PFKFB3', 'PFKFB4', 'ALDOA', 'ALDOB', 'ALDOC', 'TPI1', 'GAPDH', 'PGK1', 'PGK2', 'PGAM1', 'PGAM2', 'PGAM5', 'ENO1', 'ENO2', 'ENO3', 'PKLR', 'PKM2', 'LDHA', 'LDHAL6B', 'LDHB', 'LDHC', 'LDHD', 'GPT1', 'GPT2', 'GCKR', 'ADPGK')))
glyc.df2 <- as.data.frame(t(glyc.df[,-1]))
colnames(glyc.df2) <- glyc.df[,1]
glyc.df2 <- rownames_to_column(glyc.df2, var = 'sample')
glyc.df2$PAM50 <- all.clinical.df2$PAM50
glyc.df2$MB1 <- all.clinical.df2$MB1.forkscale.fork_0.8


# combine categories and calculate mean by each category for each enzyme expression
glyc.df2.means <- glyc.df2 %>%unite(combined, PAM50, MB1, remove = FALSE) %>% 
  pivot_longer(c(PAM50, MB1, combined), names_to = "cat_key", values_to = "cat") %>% 
  group_by(cat) %>% 
  summarise_at(vars(ADPGK:TPI1), mean, na.rm = TRUE) %>% 
  column_to_rownames(var = 'cat')

cv.metab <- data.frame(df.list$MB1.all.metab.genes.cv, df.list$MB2.all.metab.genes.cv[,2], df.list$MB3.all.metab.genes.cv[,2], stringsAsFactors = FALSE)
colnames(cv.metab) <- c('Gene', 'MB1_CV', 'MB2_CV', 'MB3_CV')

glyc.cv <- subset(cv.metab, Gene %in% c('GCK', 'HK1', 'HK2', 'HK3', 'GPI1', 'PFKL', 'PFKM', 'PFKP', 'PFKFB1', 'PFKFB2', 'PFKFB3', 'PFKFB4', 'ALDOA', 'ALDOB', 'ALDOC', 'TPI1', 'GAPDH', 'PGK1', 'PGK2', 'PGAM1', 'PGAM2', 'PGAM5', 'ENO1', 'ENO2', 'ENO3', 'PKLR', 'PKM2', 'LDHA', 'LDHAL6B', 'LDHB', 'LDHC', 'LDHD', 'GPT1', 'GPT2', 'GCKR', 'ADPGK'))
glyc.cv.MB1 <- glyc.cv[order(-glyc.cv$MB1_CV),]

glyc.df2.means <- glyc.df2.means[glyc.cv.MB1$Gene]

#  for heatmap: compare MB1_UF, MB1_LF with Pam50 categories

glyc.df2.means.heatmap <- glyc.df2.means[c('MB1_UF', 'MB1_LF', "LumA", "LumA_MB1_LF", "LumA_MB1_none", "LumB", "LumB_MB1_UF", "LumB_MB1_none"),]


# ComplexHeatmap with CVs

row_ha = rowAnnotation("MB1\nCVs" = anno_barplot(glyc.cv.MB1$MB1_CV, border = T, bar_width = 0.5, annotation_name_gp = gpar(fontface = "bold")))

Heatmap(t(glyc.df2.means.heatmap),
        name = "METABRIC \nnormalised \ngene expression",
        row_order = colnames(glyc.df2.means.heatmap),
        row_names_side = "left",
        width = unit(4, "cm"), height = unit(12, "cm"),
        column_km = 3, 
        column_title = c("UF", "PAM50", "LF"),
        right_annotation = row_ha,
        show_row_names = TRUE
)

# including all groups
glyc.df2.allgroups <- glyc.df2

glyc.df2.allgroups$MB2 <- all.clinical.df2$MB2.forkscale.fork_0.8
glyc.df2.allgroups$MB3 <- all.clinical.df2$MB3.forkscale.fork_0.8

glyc.df2.means.allgroups <- glyc.df2.allgroups %>%unite(combined, PAM50, MB1, MB2, MB3, sep = " - ", remove = FALSE) %>% 
  pivot_longer(c(PAM50, MB1, MB2, MB3, combined), names_to = "cat_key", values_to = "cat") %>% 
  group_by(cat) %>% 
  summarise_at(vars(ADPGK:TPI1), mean, na.rm = TRUE) %>% 
  column_to_rownames(var = 'cat')

# make the dataframe only with all combinations and without the NC group

glyc.df2.means.allgroups.allcombinations <- glyc.df2.means.allgroups %>% 
  rownames_to_column(var = "cat") %>% 
  filter(nchar(as.character(cat)) > 10) %>% 
  filter(!str_detect(cat, "NC"))

# create separate columns for PAM50 and the three biclusters for annotation of the heatmap


glyc.df2.means.allgroups.allcombinations.annotations <- glyc.df2.means.allgroups.allcombinations %>% 
  mutate('PAM50' = ifelse(grepl('Basal', `cat`), 'Basal',
                          ifelse(grepl('LumA', `cat`), 'LumA',
                                 ifelse(grepl('LumB', `cat`), 'LumB',
                                        ifelse(grepl('Her2', `cat`), 'Her2',
                                               ifelse(grepl('Normal', `cat`), 'Normal', 'NC')))))) %>% 
  mutate('MB1' = ifelse(grepl('MB1_UF', `cat`), 'UF',
                        ifelse(grepl('MB1_LF', `cat`), 'LF',
                               ifelse(grepl('MB1_none', `cat`), 'none','NC')))) %>% 
  mutate('MB2' = ifelse(grepl('MB2_UF', `cat`), 'UF',
                        ifelse(grepl('MB2_LF', `cat`), 'LF',
                               ifelse(grepl('MB2_none', `cat`), 'none','NC')))) %>% 
  mutate('MB3' = ifelse(grepl('MB3_UF', `cat`), 'UF',
                        ifelse(grepl('MB3_LF', `cat`), 'LF',
                               ifelse(grepl('MB3_none', `cat`), 'none','NC')))) %>% 
  column_to_rownames(var = 'cat')

# glyc.annotations <- glyc.df2.means.allgroups.allcombinations.annotations[,c("PAM50","MB1","MB2","MB3")]


fork_col <- structure(names = c("UF", "LF", 'none'), c('#bd0026','#253494', "ghostwhite"))
PAM50.col <- structure(names = c("Basal", "LumA", "LumB", "Her2", "Normal"), brewer.pal(5, "Set1"))


glyc.main <- as.matrix(glyc.df2.means.allgroups.allcombinations %>% column_to_rownames(var = 'cat'))

# glyc heatmaps ----
# the following heatmap can be plotted either with row clustering (this revelas grouping by biclusters) or by the default rwo order
# by Pam50 categories (this shows that Pam50 gives no information on metabolic patterns)

ht_list.glyc = Heatmap(glyc.main,
                      name = "METABRIC \nnormalised \ngene expression",
                      width = unit(10, "cm"), height = unit(12, "cm"),
                      # row_order = glyc.df2.means.allgroups.allcombinations$cat,
                      column_names_gp = gpar(fontsize = 6),
                      clustering_distance_columns = "pearson",
                      clustering_distance_rows = "pearson",
                      clustering_method_rows = "ward.D2",
                      clustering_method_columns = "ward.D2",
                      row_names_side = "left",
                      column_names_side = "top",
                      show_row_names = FALSE,
                      show_column_names = TRUE,
                      row_title = NULL,
                      column_title = NULL,
                      column_split = 3,
                      # column_order = colnames(glyc.df2.means.heatmap)
                      row_km = 4
) +
  Heatmap(glyc.df2.means.allgroups.allcombinations.annotations["PAM50"],
          name = "PAM50",
          width = unit(0.3, "cm"), height = unit(20, "cm"),
          col = PAM50.col,
          show_row_names = FALSE,
          column_names_gp = gpar(fontsize = 9),
          column_names_side = "top") +
  Heatmap(glyc.df2.means.allgroups.allcombinations.annotations[,c("MB1","MB2","MB3")],
          name = "biclusters\nforks",
          width = unit(0.9, "cm"), height = unit(20, "cm"),
          col = fork_col,
          column_names_gp = gpar(fontsize = 9),
          show_row_names = FALSE,
          column_names_side = "top",
          cluster_columns = FALSE)

pdf("...", width = 8, height = 10)

draw(ht_list.glyc)

dev.off()

ht_glyc <- Heatmap(glyc.main,
                  name = "METABRIC \nnormalised \ngene expression",
                  width = unit(10, "cm"), height = unit(20, "cm"),
                  # row_order = glyc.df2.means.allgroups.allcombinations$cat,
                  clustering_distance_columns = "pearson",
                  clustering_distance_rows = "pearson",
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2",
                  row_names_side = "left",
                  column_names_side = "top",
                  show_row_names = FALSE,
                  show_column_names = FALSE,
                  row_title = NULL,
                  column_title = NULL,
                  column_split = 4,
                  # column_order = colnames(glyc.df2.means.heatmap)
                  row_km = 3)

column_order(ht_glyc)


glyc.clusters.fork.list <- list()
glyc.clusters.fork.list$"slice1" <- colnames(glyc.main)[column_order(ht_glyc)[[1]]]
glyc.clusters.fork.list$"slice2" <- colnames(glyc.main)[column_order(ht_glyc)[[2]]]
glyc.clusters.fork.list$"slice3" <- colnames(glyc.main)[column_order(ht_glyc)[[3]]]
glyc.clusters.fork.list$"slice4" <- colnames(glyc.main)[column_order(ht_glyc)[[4]]]

