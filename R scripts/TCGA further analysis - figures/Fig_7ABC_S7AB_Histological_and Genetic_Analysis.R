
library(readr)
library(dplyr)
library(tidyverse)
library(MASS)
library(MCbiclust)
library(ggplot2)
library(ComplexHeatmap)
library(ggthemes)


# 1. Load TCGA RNASeq data - assemble all biclusters ----

load('../TCGA_MB1_RNASeq_data_nonorm.RData')
load('../TCGA_MB2_RNASeq_data_nonorm.RData')
load('../TCGA_MB3_RNASeq_data_nonorm.RData')


TCGA.all.biclusters.RNAseq.df <- TCGA.MB1.RNAseq.df %>% left_join(TCGA.MB2.RNAseq.df) %>% left_join(TCGA.MB3.RNAseq.df)
TCGA.all.biclusters.RNAseq.df <- TCGA.all.biclusters.RNAseq.df[,c(1:3, 7:10, 4:6)]

# 2. Load and combine new clinical and histological data with bicluster data ----
consensus <- read_csv("consensus.csv")
TCGA_molecular_profiles <- read_csv("TCGA_molecular_profiles.csv")

TCGA.all.biclusters.RNAseq.df$sample <- substr(TCGA.all.biclusters.RNAseq.df$X,1,12)
TCGA.all.biclusters.RNAseq.df <- TCGA.all.biclusters.RNAseq.df[,c(1,11,2:10)]
TCGA.all.biclusters.RNAseq.df <- na.omit(TCGA.all.biclusters.RNAseq.df)
TCGA.all.biclusters.RNAseq.df$MB1.fork <- ForkClassifier(TCGA.all.biclusters.RNAseq.df$MB1.pc1,300)
TCGA.all.biclusters.RNAseq.df <- TCGA.all.biclusters.RNAseq.df[,c(1:4,12,5:11)]
TCGA.all.biclusters.RNAseq.df <- TCGA.all.biclusters.RNAseq.df[order(TCGA.all.biclusters.RNAseq.df$MB3.index),]
TCGA.all.biclusters.RNAseq.df$MB3.fork <- ForkClassifier(TCGA.all.biclusters.RNAseq.df$MB3.pc1,300)
TCGA.all.biclusters.RNAseq.df <- TCGA.all.biclusters.RNAseq.df[,c(1:7,13,8:12)]
TCGA.all.biclusters.RNAseq.df <- TCGA.all.biclusters.RNAseq.df[order(TCGA.all.biclusters.RNAseq.df$MB2.index),]
TCGA.all.biclusters.RNAseq.df$MB2.fork <- ForkClassifier(TCGA.all.biclusters.RNAseq.df$MB2.pc1,300)
TCGA.all.biclusters.RNAseq.df <- TCGA.all.biclusters.RNAseq.df[,c(1:10,14,11:13)]

TCGA_molecular_profiles$sample <- gsub('\\.','-',TCGA_molecular_profiles$Row.names)

TCGA.all.biclusters.RNAseq.df2 <- TCGA.all.biclusters.RNAseq.df %>%
  inner_join(consensus,by = 'sample') %>%
  inner_join(TCGA_molecular_profiles,by = 'sample')

TCGA.all.biclusters.RNAseq.df2$Tumor.Purity <- as.numeric(TCGA.all.biclusters.RNAseq.df2$Tumor.Purity)
TCGA.all.biclusters.RNAseq.df2$MB1.fork <- factor(TCGA.all.biclusters.RNAseq.df2$MB1.fork,
                                    levels = c('None','Lower','Upper'))
TCGA.all.biclusters.RNAseq.df2$MB3.fork <- factor(TCGA.all.biclusters.RNAseq.df2$MB3.fork,
                                                   levels = c('None','Lower','Upper'))
TCGA.all.biclusters.RNAseq.df2$MB2.fork <- factor(TCGA.all.biclusters.RNAseq.df2$MB2.fork,
                                                   levels = c('None','Lower','Upper'))

# add thresholdforks
TCGA.all.biclusters.RNAseq.df2 <- mutate(TCGA.all.biclusters.RNAseq.df2, MB1.forkscale = MB1.pc1/MB1.index)
TCGA.all.biclusters.RNAseq.df2 <- mutate(TCGA.all.biclusters.RNAseq.df2, MB3.forkscale = MB3.pc1/MB3.index)
TCGA.all.biclusters.RNAseq.df2 <- mutate(TCGA.all.biclusters.RNAseq.df2, MB2.forkscale = MB2.pc1/MB2.index)

TCGA.all.biclusters.RNAseq.df2 <- mutate(TCGA.all.biclusters.RNAseq.df2, MB1.forkscale.fork_0.4 = ifelse(`MB1.forkscale` <= -0.04 , "Lower", ifelse(`MB1.forkscale` >= 0.04 , "Upper", "none")))
TCGA.all.biclusters.RNAseq.df2 <- mutate(TCGA.all.biclusters.RNAseq.df2, MB1.forkscale.fork_0.8 = ifelse(`MB1.forkscale` <= -0.08 , "Lower", ifelse(`MB1.forkscale` >= 0.08 , "Upper", "none")))
TCGA.all.biclusters.RNAseq.df2 <- mutate(TCGA.all.biclusters.RNAseq.df2, MB3.forkscale.fork_0.4 = ifelse(`MB3.forkscale` <= -0.04 , "Lower", ifelse(`MB3.forkscale` >= 0.04 , "Upper", "none")))
TCGA.all.biclusters.RNAseq.df2 <- mutate(TCGA.all.biclusters.RNAseq.df2, MB3.forkscale.fork_0.8 = ifelse(`MB3.forkscale` <= -0.08 , "Lower", ifelse(`MB3.forkscale` >= 0.08 , "Upper", "none")))
TCGA.all.biclusters.RNAseq.df2 <- mutate(TCGA.all.biclusters.RNAseq.df2, MB2.forkscale.fork_0.4 = ifelse(`MB2.forkscale` <= -0.04 , "Lower", ifelse(`MB2.forkscale` >= 0.04 , "Upper", "none")))
TCGA.all.biclusters.RNAseq.df2 <- mutate(TCGA.all.biclusters.RNAseq.df2, MB2.forkscale.fork_0.8 = ifelse(`MB2.forkscale` <= -0.08 , "Lower", ifelse(`MB2.forkscale` >= 0.08 , "Upper", "none")))


TCGA.all.biclusters.RNAseq.df2 <- mutate(TCGA.all.biclusters.RNAseq.df2, `Hist_type` = ifelse(grepl("Ductal", `hist_type`), "IDC",
                                                                ifelse(grepl("Lobular", `hist_type`), "ILC",
                                                                       ifelse(grepl("Mixed", `hist_type`), "mixed", "other"))))

TCGA.all.biclusters.RNAseq.df2 <- mutate(TCGA.all.biclusters.RNAseq.df2, `Nuclear_pleo` = ifelse(grepl("Small", `nuc_pleo`), "regular",
                                                                                              ifelse(grepl("Increase", `nuc_pleo`), "moderate",
                                                                                                     ifelse(grepl("Marked", `nuc_pleo`), "marked", "other"))))

TCGA.all.biclusters.RNAseq.df2 <- mutate(TCGA.all.biclusters.RNAseq.df2, `mitosisS` = ifelse(grepl("Low", `mitosis`), "low",
                                                                                                 ifelse(grepl("Medium", `mitosis`), "medium",
                                                                                                        ifelse(grepl("High", `mitosis`), "high", "other"))))

TCGA.all.biclusters.RNAseq.df2 <- mutate(TCGA.all.biclusters.RNAseq.df2, `epi_areaS` = ifelse(grepl("Low", `epi_area`), "low",
                                                                                             ifelse(grepl("Moderate", `epi_area`), "moderate",
                                                                                                    ifelse(grepl("High", `epi_area`), "high", "other"))))




TCGA.all.biclusters.RNAseq.df2 <- TCGA.all.biclusters.RNAseq.df2 %>% 
  mutate(MB1.forkscale.log = MB1.pc1/(log(MB1.index))) %>% 
  mutate(MB2.forkscale.log = MB2.pc1/(log(MB2.index))) %>% 
  mutate(MB3.forkscale.log = MB3.pc1/(log(MB3.index)))


save(TCGA.all.biclusters.RNAseq.df2, file='TCGA.all.biclusters.RNAseq.Rdata')

WriteXLS(TCGA.all.biclusters.RNAseq.df2, ExcelFileName = 'TCGA-all-biclusters-RNAseq.xlsx')


# 3. Create 'three-way' plots Fig 7B, S7A: MB1 vs MB2 ----

theme_fork3 <- theme_tufte(base_size = 16, base_family = 'sans') +
  theme(aspect.ratio = 1/1) +
  theme(axis.ticks.length = unit(8, "pt")) +
  theme(axis.ticks = element_line(colour = "black", size = 0.6)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 2, r = 0, b = 0, l = 0))) +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.line = element_blank())+
  theme(axis.text.y = element_text(color = "black", size = 16))+
  theme(axis.text.x = element_text(color = "black", size = 16))+
  theme(legend.margin=margin(0,0,0,0), legend.box.margin=margin(20,0,0,-310), legend.justification="top")+
  theme(panel.grid.minor = element_blank())+
  theme(legend.title=element_text(size=14))+
  theme(legend.title.align = 0.5) +
  theme(plot.title = element_text(size = 15, face = 'plain', hjust = 0.1, vjust = -20))

theme_set(theme_fork3)
grids <- data.frame(c(0,0), c(-10,10))
colnames(grids) <- c('x','y')

# hist type

pdf("..", height = 5.08, width = 10)

ggplot(TCGA.all.biclusters.RNAseq.df2, aes(MB1.forkscale.log, MB2.forkscale.log)) +
  # geom_point(data = TCGA.all.biclusters.RNAseq.df2, aes(MB1.forkscale, MB2.forkscale), colour = '#d9d9d9', alpha=0.1, shape = 16) +
  geom_point(aes(colour = Hist_type), alpha=0.7, size=2, shape = 16) +
  scale_colour_manual(values = c("ILC" = "#e41a1c", "mixed"="#fec44f", "IDC"="#4daf4a", "other"="#bdbdbd"), name = 'histological type') +
  geom_rangeframe(data=data.frame(MB1.forkscale.log=c(-10,10), MB2.forkscale.log=c(-10,10)), size = 1.2) +
  scale_x_continuous(breaks = seq(-10, 10, 5), limits = c(-10,10), minor_breaks = seq(-20,20,1)) +
  scale_y_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 0.8))) +
  geom_line(data = grids, aes(x=x, y=y), linetype = "dotted") +
  geom_line(data = grids, aes(x=y, y=x), linetype = "dotted") +
  labs(x ="MB1 forkscale", y = "MB2 forkscale")

dev.off()


# Nuclear_pleo

pdf("..", height = 5.08, width = 10)

ggplot(TCGA.all.biclusters.RNAseq.df2, aes(MB1.forkscale.log, MB2.forkscale.log)) +
  # geom_point(data = TCGA.all.biclusters.RNAseq.df2, aes(MB1.forkscale, MB2.forkscale), colour = '#d9d9d9', alpha=0.1, shape = 16) +
  geom_point(aes(colour = Nuclear_pleo), alpha=0.7, size=2, shape = 16) +
  scale_colour_manual(values = c("marked" = "#e41a1c", "moderate"="#4daf4a", "small"="#fec44f"), name = 'nuclear pleomorphism') +
  geom_rangeframe(data=data.frame(MB1.forkscale.log=c(-10,10), MB2.forkscale.log=c(-10,10)), size = 1.2) +
  scale_x_continuous(breaks = seq(-10, 10, 5), limits = c(-10,10), minor_breaks = seq(-20,20,1)) +
  scale_y_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 0.8))) +
  geom_line(data = grids, aes(x=x, y=y), linetype = "dotted") +
  geom_line(data = grids, aes(x=y, y=x), linetype = "dotted") +
  labs(x ="MB1 forkscale", y = "MB2 forkscale")

dev.off()

# mitosis

pdf("..", height = 5.08, width = 10)

ggplot(TCGA.all.biclusters.RNAseq.df2, aes(MB1.forkscale.log, MB2.forkscale.log)) +
  # geom_point(data = TCGA.all.biclusters.RNAseq.df2, aes(MB1.forkscale, MB2.forkscale), colour = '#d9d9d9', alpha=0.1, shape = 16) +
  geom_point(aes(colour = mitosisS), alpha=0.7, size=2, shape = 16) +
  scale_colour_manual(values = c("high" = "#e41a1c", "medium"="#4daf4a", "low"="#fec44f", "other"="#bdbdbd"), name = 'mitosis') +
  geom_rangeframe(data=data.frame(MB1.forkscale.log=c(-10,10), MB2.forkscale.log=c(-10,10)), size = 1.2) +
  scale_x_continuous(breaks = seq(-10, 10, 5), limits = c(-10,10), minor_breaks = seq(-20,20,1)) +
  scale_y_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 0.8))) +
  geom_line(data = grids, aes(x=x, y=y), linetype = "dotted") +
  geom_line(data = grids, aes(x=y, y=x), linetype = "dotted") +
  labs(x ="MB1 forkscale", y = "MB2 forkscale")

dev.off()


# epi_area

pdf("..", height = 5.08, width = 10)

ggplot(TCGA.all.biclusters.RNAseq.df2, aes(MB1.forkscale.log, MB2.forkscale.log)) +
  # geom_point(data = TCGA.all.biclusters.RNAseq.df2, aes(MB1.forkscale, MB2.forkscale), colour = '#d9d9d9', alpha=0.1, shape = 16) +
  geom_point(aes(colour = epi_areaS), alpha=0.7, size=2, shape = 16) +
  scale_colour_manual(values = c("high" = "#e41a1c", "moderate"="#4daf4a", "low"="#fec44f", "other"="#bdbdbd"), name = 'epi area') +
  geom_rangeframe(data=data.frame(MB1.forkscale.log=c(-10,10), MB2.forkscale.log=c(-10,10)), size = 1.2) +
  scale_x_continuous(breaks = seq(-10, 10, 5), limits = c(-10,10), minor_breaks = seq(-20,20,1)) +
  scale_y_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 0.8))) +
  geom_line(data = grids, aes(x=x, y=y), linetype = "dotted") +
  geom_line(data = grids, aes(x=y, y=x), linetype = "dotted") +
  labs(x ="MB1 forkscale", y = "MB2 forkscale")

dev.off()


# inflam

pdf("..", height = 5.08, width = 10)

ggplot(TCGA.all.biclusters.RNAseq.df2, aes(MB1.forkscale.log, MB2.forkscale.log)) +
  # geom_point(data = TCGA.all.biclusters.RNAseq.df2, aes(MB1.forkscale, MB2.forkscale), colour = '#d9d9d9', alpha=0.1, shape = 16) +
  geom_point(aes(colour = inflam), alpha=0.7, size=2, shape = 16) +
  scale_colour_manual(values = c("Present" = "#e41a1c", "Absent"="#fec44f", "NA"="#bdbdbd"), name = 'inflammation') +
  geom_rangeframe(data=data.frame(MB1.forkscale.log=c(-10,10), MB2.forkscale.log=c(-10,10)), size = 1.2) +
  scale_x_continuous(breaks = seq(-10, 10, 5), limits = c(-10,10), minor_breaks = seq(-20,20,1)) +
  scale_y_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 0.8))) +
  geom_line(data = grids, aes(x=x, y=y), linetype = "dotted") +
  geom_line(data = grids, aes(x=y, y=x), linetype = "dotted") +
  labs(x ="MB1 forkscale", y = "MB2 forkscale")

dev.off()



# 4. plot functions for MB1, to generate Fig 7A,C ---- 

hist.plot.fun.MB1 <- function(x){
  p <- ggplot(TCGA.all.biclusters.RNAseq.df2 , aes(MB1.index, MB1.pc1)) +
    geom_point(aes(colour = eval(parse(text = x)))) + 
    scale_colour_discrete(name = x)
  
  print(p)
  tbl = table(TCGA.all.biclusters.RNAseq.df2$MB1.fork, TCGA.all.biclusters.RNAseq.df2[[x]]) 
  print(tbl)
  print(chisq.test(tbl))
}

hist.plot.fun2.MB1 <- function(x){
  p <- ggplot(TCGA.all.biclusters.RNAseq.df2 , aes(MB1.index, MB1.pc1)) +
    geom_point(aes(colour = eval(parse(text = x)))) + 
    scale_colour_continuous(name = x)
  p1 <- ggplot(TCGA.all.biclusters.RNAseq.df2, aes(MB1.fork,  eval(parse(text = x)))) + 
    geom_boxplot() + 
    ylab(x)
  
  print(p)
  print(p1)
  
  fit <- aov(eval(parse(text = x)) ~ MB1.fork, data=TCGA.all.biclusters.RNAseq.df2)
  print(summary(fit))
  print(TukeyHSD(fit))
}

# Two different outputs for PAM50
# Calculated from TCGA RNASeq data directly:
hist.plot.fun.MB1('PAM50.x')
# From histological TCGA data set:
hist.plot.fun.MB1('PAM50.y')


# discrete variables:
hist.plot.fun.MB1('epi_area')
hist.plot.fun.MB1('inflam')
hist.plot.fun.MB1('lcis')
hist.plot.fun.MB1('apo_feat')
hist.plot.fun.MB1('dcis')
hist.plot.fun.MB1('epi_tube')
hist.plot.fun.MB1('lymp')
hist.plot.fun.MB1('necrosis')
hist.plot.fun.MB1('nuc_pleo')
hist.plot.fun.MB1('fib_focus')
hist.plot.fun.MB1('mitosis')
hist.plot.fun.MB1('hist_type')
hist.plot.fun.MB1('miRNA.Clusters')
hist.plot.fun.MB1('methylation.Clusters')
hist.plot.fun.MB1('RPPA.Clusters')
hist.plot.fun.MB1('CN.Clusters')

# continuous variables:
hist.plot.fun2.MB1('Proliferation')
hist.plot.fun2.MB1('Tumor.Purity')

# 5. plot functions for MB2, to generate Fig 7A,C ---- 

hist.plot.fun.MB2 <- function(x){
  p <- ggplot(TCGA.all.biclusters.RNAseq.df2 , aes(MB2.index, MB2.pc1)) +
    geom_point(aes(colour = eval(parse(text = x)))) + 
    scale_colour_discrete(name = x)
  
  print(p)
  tbl = table(TCGA.all.biclusters.RNAseq.df2$MB2.fork, TCGA.all.biclusters.RNAseq.df2[[x]]) 
  print(tbl)
  print(chisq.test(tbl))
}

hist.plot.fun2 <- function(x){
  p <- ggplot(TCGA.all.biclusters.RNAseq.df2 , aes(MB2.index, MB2.pc1)) +
    geom_point(aes(colour = eval(parse(text = x)))) + 
    scale_colour_continuous(name = x)
  p1 <- ggplot(TCGA.all.biclusters.RNAseq.df2, aes(MB2.fork,  eval(parse(text = x)))) + 
    geom_boxplot() + 
    ylab(x)
  
  print(p)
  print(p1)
  
  fit <- aov(eval(parse(text = x)) ~ MB2.fork, data=TCGA.all.biclusters.RNAseq.df2)
  print(summary(fit))
  print(TukeyHSD(fit))
}

# Two different outputs for PAM50
# Calculated from TCGA RNASeq data directly:
hist.plot.fun.MB2('PAM50.x')
# From histological TCGA data set:
hist.plot.fun.MB2('PAM50.y')

# discrete variables:
hist.plot.fun.MB2('epi_area')
hist.plot.fun.MB2('inflam')
hist.plot.fun.MB2('lcis')
hist.plot.fun.MB2('apo_feat')
hist.plot.fun.MB2('dcis')
hist.plot.fun.MB2('epi_tube')
hist.plot.fun.MB2('lymp')
hist.plot.fun.MB2('necrosis')
hist.plot.fun.MB2('nuc_pleo')
hist.plot.fun.MB2('fib_focus')
hist.plot.fun.MB2('mitosis')
hist.plot.fun.MB2('hist_type')
hist.plot.fun.MB2('miRNA.Clusters')
hist.plot.fun.MB2('methylation.Clusters')
hist.plot.fun.MB2('RPPA.Clusters')
hist.plot.fun.MB2('CN.Clusters')

# continuous variables:
hist.plot.fun2('Proliferation')
hist.plot.fun2('Tumor.Purity')

# 6. plot functions for MB3, to generate Fig S7B ----

hist.plot.fun.MB3 <- function(x){
  p <- ggplot(TCGA.all.biclusters.RNAseq.df2 , aes(MB3.index, MB3.pc1)) +
    geom_point(aes(colour = eval(parse(text = x)))) + 
    scale_colour_discrete(name = x)
  
  print(p)
  tbl = table(TCGA.all.biclusters.RNAseq.df2$MB3.fork, TCGA.all.biclusters.RNAseq.df2[[x]]) 
  print(tbl)
  print(chisq.test(tbl))
}

hist.plot.fun2 <- function(x){
  p <- ggplot(TCGA.all.biclusters.RNAseq.df2 , aes(MB3.index, MB3.pc1)) +
    geom_point(aes(colour = eval(parse(text = x)))) + 
    scale_colour_continuous(name = x)
  p1 <- ggplot(TCGA.all.biclusters.RNAseq.df2, aes(MB3.fork,  eval(parse(text = x)))) + 
    geom_boxplot() + 
    ylab(x)
  
  print(p)
  print(p1)
  
  fit <- aov(eval(parse(text = x)) ~ MB3.fork, data=TCGA.all.biclusters.RNAseq.df2)
  print(summary(fit))
  print(TukeyHSD(fit))
}

# Two different outputs for PAM50
# Calculated from TCGA RNASeq data directly:
hist.plot.fun.MB3('PAM50.x')
# From histological TCGA data set:
hist.plot.fun.MB3('PAM50.y')

# discrete variables:
hist.plot.fun.MB3('epi_area')
hist.plot.fun.MB3('inflam')
hist.plot.fun.MB3('lcis')
hist.plot.fun.MB3('apo_feat')
hist.plot.fun.MB3('dcis')
hist.plot.fun.MB3('epi_tube')
hist.plot.fun.MB3('lymp')
hist.plot.fun.MB3('necrosis')
hist.plot.fun.MB3('nuc_pleo')
hist.plot.fun.MB3('fib_focus')
hist.plot.fun.MB3('mitosis')
hist.plot.fun.MB3('hist_type')
hist.plot.fun.MB3('miRNA.Clusters')
hist.plot.fun.MB3('methylation.Clusters')
hist.plot.fun.MB3('RPPA.Clusters')
hist.plot.fun.MB3('CN.Clusters')

# continuous variables:
hist.plot.fun2('Proliferation')
hist.plot.fun2('Tumor.Purity')


# 7. heatmap Fig. 7A chi-squared p values ----

qsmatrix <- read_csv('qsmatrix.csv')
qsmatrix <- column_to_rownames(qsmatrix, var = '...1')
qsmatrix <- as.matrix(qsmatrix)

col.pal <- RColorBrewer::brewer.pal(9, "Reds")

robust_dist = function(x, y) {
  qx = quantile(x, c(0.1, 0.9))
  qy = quantile(y, c(0.1, 0.9))
  l = x > qx[1] & x < qx[2] & y > qy[1] & y < qy[2]
  x = x[l]
  y = y[l]
  sqrt(sum((x - y)^2))
}

pdf("..", onefile = TRUE, width = 5, height = 14)

Heatmap(qsmatrix,
        name = '-logP',
        col = col.pal,
        clustering_distance_rows = robust_dist,
        clustering_distance_columns = robust_dist,
        # clustering_method_rows = 'ward.D2',
        row_dend_reorder = FALSE,
        column_order = order(as.numeric(gsub("MB", "", colnames(qsmatrix)))),
        row_km = 3,
        column_names_side = "top",
        column_names_rot = 0,
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 8),
        # row_title = "cluster_%s",
        # column_km = 2,
        # column_split = 3,
        width = unit(3, "cm"), height = unit(9, "cm")
)

dev.off()



# 8. Heatmaps of prop tables Fig 7C, S7B ----

# prop tables ----

discrete.variables <- c('PAM50.y','epi_area', 'inflam', 'lcis', 'apo_feat', 'dcis', 'epi_tube','lymp',
                        'necrosis', 'nuc_pleo', 'fib_focus', 'mitosis', 'hist_type', 'miRNA.Clusters',
                        'methylation.Clusters', 'RPPA.Clusters', 'CN.Clusters')

# MB1
prop.table.list <- list()
for(i in discrete.variables){
  prop.table.list[[i]] <- prop.table(table(TCGA.all.biclusters.RNAseq.df2$MB1.fork, TCGA.all.biclusters.RNAseq.df2[[i]]),1)
  colnames(prop.table.list[[i]]) <- paste(i, colnames(prop.table.list[[i]]), sep = ' : ')
}
prop.table.all.MB1 <- Reduce(cbind,prop.table.list)
prop.table.all.minus.MB1 <-prop.table.all.MB1[,-grep('Absent',colnames(prop.table.all.MB1))]
prop.table.all.minus.MB1 <-t(prop.table.all.minus.MB1)[-grep('Missing',colnames(prop.table.all.minus.MB1)),]

# MB2
prop.table.list <- list()
for(i in discrete.variables){
  prop.table.list[[i]] <- prop.table(table(TCGA.all.biclusters.RNAseq.df2$MB2.fork, TCGA.all.biclusters.RNAseq.df2[[i]]),1)
  colnames(prop.table.list[[i]]) <- paste(i, colnames(prop.table.list[[i]]), sep = ' : ')
}
prop.table.all.MB2 <- Reduce(cbind,prop.table.list)
prop.table.all.minus.MB2 <-prop.table.all.MB2[,-grep('Absent',colnames(prop.table.all.MB2))]
prop.table.all.minus.MB2 <-t(prop.table.all.minus.MB2)[-grep('Missing',colnames(prop.table.all.minus.MB2)),]

# MB3
prop.table.list <- list()
for(i in discrete.variables){
  prop.table.list[[i]] <- prop.table(table(TCGA.all.biclusters.RNAseq.df2$MB3.fork, TCGA.all.biclusters.RNAseq.df2[[i]]),1)
  colnames(prop.table.list[[i]]) <- paste(i, colnames(prop.table.list[[i]]), sep = ' : ')
}
prop.table.all.MB3 <- Reduce(cbind,prop.table.list)
prop.table.all.minus.MB3 <-prop.table.all.MB3[,-grep('Absent',colnames(prop.table.all.MB3))]
prop.table.all.minus.MB3 <-t(prop.table.all.minus.MB3)[-grep('Missing',colnames(prop.table.all.minus.MB3)),]

# all forks

prop.table.all <- rownames_to_column(as.data.frame(prop.table.all.minus.MB1), 'rn') %>% 
  left_join(rownames_to_column(as.data.frame(prop.table.all.minus.MB2),'rn'), by = 'rn') %>% 
  left_join(rownames_to_column(as.data.frame(prop.table.all.minus.MB3),'rn'), by = 'rn')


prop.table.all <- prop.table.all[,c(1,3,4,6,7,9,10)]
colnames(prop.table.all) <- c('rn','MB1_LF', 'MB1_UF', 'MB2_LF', 'MB2_UF', 'MB3_LF', 'MB3_UF')
prop.table.all.df <- prop.table.all[,-1]
rownames(prop.table.all.df) <- prop.table.all[,1]
prop.table.all.matrix <- as.matrix(prop.table.all.df)
prop.table.MB1_2.matrix <- prop.table.all.matrix[,c(1:4)]

# heatmap MB1 MB2 Fig 7C ----

col.pal <- RColorBrewer::brewer.pal(9, "Reds")

robust_dist = function(x, y) {
  qx = quantile(x, c(0.1, 0.9))
  qy = quantile(y, c(0.1, 0.9))
  l = x > qx[1] & x < qx[2] & y > qy[1] & y < qy[2]
  x = x[l]
  y = y[l]
  sqrt(sum((x - y)^2))
}

pdf("..", onefile = TRUE, width = 10, height = 14)

Heatmap(prop.table.all.matrix[,c(1:4)],
        name = 'fraction \nof samples',
        col = col.pal,
        clustering_distance_rows = robust_dist,
        clustering_distance_columns = robust_dist,
        # clustering_method_rows = 'ward.D2',
        row_dend_reorder = TRUE,
        row_km = 6,
        row_title = "cluster_%s",
        show_row_names = TRUE,
        # column_km = 2,
        # column_split = 3,
        width = unit(2.5, "cm"), height = unit(14, "cm")
)

# heatmap MB3 Fig S7B ----

col.pal <- RColorBrewer::brewer.pal(9, "Reds")

robust_dist = function(x, y) {
  qx = quantile(x, c(0.1, 0.9))
  qy = quantile(y, c(0.1, 0.9))
  l = x > qx[1] & x < qx[2] & y > qy[1] & y < qy[2]
  x = x[l]
  y = y[l]
  sqrt(sum((x - y)^2))
}

pdf("..", onefile = TRUE, width = 10, height = 14)

Heatmap(prop.table.all.matrix[,c(5,6)],
        name = 'fraction \nof samples',
        col = col.pal,
        clustering_distance_rows = robust_dist,
        clustering_distance_columns = robust_dist,
        # clustering_method_rows = 'ward.D2',
        row_dend_reorder = TRUE,
        row_km = 6,
        row_title = "cluster_%s",
        show_row_names = TRUE,
        # column_km = 2,
        # column_split = 3,
        width = unit(2.5, "cm"), height = unit(14, "cm")
        
)


dev.off()

