
library(gdata)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(ggnewscale)
library(readxl)


# Load data
load("../METABRIC_starting_data_final_corr_groups.Rdata")

metab.genes.xls <- read_excel('GS metabolic genes list.xlsx', sheet=1)

cv.diff.metab <- data.frame(df.list$MB1.all.metab.genes.cv, df.list$MB2.all.metab.genes.cv[,2], ((df.list$MB1.all.metab.genes.cv[,2]) - (df.list$MB2.all.metab.genes.cv[,2])), stringsAsFactors = FALSE)
colnames(cv.diff.metab) <- c('Gene', 'MB1_CV', 'MB2_CV', 'MB1-MB2')

# glutamine metabolism ----

glut.df <- subset(cv.diff.metab, Gene %in% c('GFPT1', 'GFPT2', 'PPAT', 'GLUL', 'GLS', 'GLS2', 'GLUD1', 'GLUD2', 'GAD1', 'GAD2', 'GCLC', 'GCLM', 'GSS', 'GGCT', 'GPT', 'GPT2', 'PSAT1'))

glut.df <- glut.df[order(glut.df$MB1_CV),]

glut.df.g <-gather(glut.df, `MB1_CV`, `MB2_CV`, key = "bicluster", value = "CV")

pdf("..", height = 4, width = 6)

ggplot(glut.df, aes(reorder(Gene, MB1_CV, decreasing = TRUE), MB1_CV, fill = MB1_CV, position = 'stack')) + 
  scale_fill_gradient2('MB1', low = '#a6bddb', mid = 'white', high = '#D55E00', midpoint = 0, aesthetics = "fill") +
  geom_col()+
  new_scale("fill")+
  # new_scale_fill() +
  geom_col(aes(reorder(Gene, MB1_CV, decreasing = TRUE), MB2_CV, fill = MB1_CV, position = 'stack')) +
  scale_fill_gradient2('MB2', low = '#980043', mid = 'white', high = '#1b7837', midpoint = 0, aesthetics = "fill")+
    theme_classic() +
  theme(axis.title.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(angle = 40, hjust = 0.7, vjust = 1.5)) +
  labs(fill = 'bicluster') +
  ylab(label = 'CV') +
  theme(legend.position="bottom", legend.title.align=0.8) +
  xlab(label = 'Genes') +
  ylim(-2, 2)

dev.off()

glut.df.g.forks <- glut.df.g %>% mutate (forks = ifelse(bicluster == 'MB1_CV' & CV > '0', 'MB1_UF', 
                                                        ifelse(bicluster == 'MB1_CV' & CV < '0', 'MB1_LF', 
                                                               ifelse(bicluster == 'MB2_CV' & CV > '0', 'MB2_UF', 'MB2_LF'))))
  
glut.df.g.forks$MB1_CV <- as.numeric(paste(glut.df$MB1_CV)) 

pdf("..", height = 4, width = 3)

ggplot(glut.df.g.forks, aes(reorder(Gene, MB1_CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB2_UF", "MB1_UF", 'MB2_LF', 'MB1_LF')), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#1b7837','#D55E00','#980043','#a6bddb'))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
        ) +
  # scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2))+
  ylab(label = 'CV value') +
  theme(legend.position="none") +
  ylim(-1.3, 2)

dev.off()


# including MB3, and create a similar plot only with that bicluster:

cv.MB3 <- data.frame(df.list$MB3.all.metab.genes.cv)

glut.cv.MB3 <- subset(cv.MB3, Gene %in% c('GFPT1', 'GFPT2', 'PPAT', 'GLUL', 'GLS', 'GLS2', 'GLUD1', 'GLUD2', 'GAD1', 'GAD2', 'GCLC', 'GCLM', 'GSS', 'GGCT', 'GPT', 'GPT2', 'PSAT1'))
glut.cv.MB3 <- glut.cv.MB3[order(glut.cv.MB3$CV),]
glut.cv.MB3 <- mutate(glut.cv.MB3, forks = ifelse(CV > 0, 'MB3_UF', 'MB3_LF'))

pdf("..", height = 4, width = 1.5)

ggplot(glut.cv.MB3, aes(reorder(Gene, CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB3_UF", "MB3_LF")), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#762a83','#d7301f'))+
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  ylab(label = 'CV value') +
  theme(legend.position="none")

dev.off()


# glycolysis ----

glyc.df <- subset(cv.diff.metab, Gene %in% c('GCK', 'HK1', 'HK2', 'HK3', 'GPI1', 'PFKL', 'PFKM', 'PFKP', 'PFKFB1', 'PFKFB2', 'PFKFB3', 'PFKFB4', 'ALDOA', 'ALDOB', 'ALDOC', 'TPI1', 'GAPDH', 'PGK1', 'PGK2', 'PGAM1', 'PGAM2', 'PGAM5', 'ENO1', 'ENO2', 'ENO3', 'PKLR', 'PKM2', 'LDHA', 'LDHAL6B', 'LDHB', 'LDHC', 'LDHD', 'GPT1', 'GPT2', 'GCKR', 'ADPGK'))
glyc.df <- glyc.df[order(glyc.df$MB1_CV),]
glyc.df.g <- gather(glyc.df, `MB1_CV`, `MB2_CV`, key = "bicluster", value = "CV")
glyc.df.g.forks <- glyc.df.g %>% mutate (forks = ifelse(bicluster == 'MB1_CV' & CV > '0', 'MB1_UF', 
                                                        ifelse(bicluster == 'MB1_CV' & CV < '0', 'MB1_LF', 
                                                               ifelse(bicluster == 'MB2_CV' & CV > '0', 'MB2_UF', 'MB2_LF'))))
glyc.df.g.forks$MB1_CV <- as.numeric(paste(glyc.df$MB1_CV)) 



pdf("..", height = 4, width = 3)


ggplot(glyc.df.g.forks, aes(reorder(Gene, MB1_CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB2_UF", "MB1_UF", 'MB2_LF', 'MB1_LF')), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#1b7837','#D55E00','#980043','#a6bddb'))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  # scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2))+
  ylab(label = 'CV value') +
  theme(legend.position="none") +
  ylim(-1.3, 2)

dev.off()


# including MB3, and create a similar plot only with that bicluster:

cv.MB3 <- data.frame(df.list$MB3.all.metab.genes.cv)

glyc.cv.MB3 <- subset(cv.MB3, Gene %in% c('GCK', 'HK1', 'HK2', 'HK3', 'GPI1', 'PFKL', 'PFKM', 'PFKP', 'PFKFB1', 'PFKFB2', 'PFKFB3', 'PFKFB4', 'ALDOA', 'ALDOB', 'ALDOC', 'TPI1', 'GAPDH', 'PGK1', 'PGK2', 'PGAM1', 'PGAM2', 'PGAM5', 'ENO1', 'ENO2', 'ENO3', 'PKLR', 'PKM2', 'LDHA', 'LDHAL6B', 'LDHB', 'LDHC', 'LDHD', 'GPT1', 'GPT2', 'GCKR', 'ADPGK'))
glyc.cv.MB3 <- glyc.cv.MB3[order(glyc.cv.MB3$CV),]
glyc.cv.MB3 <- mutate(glyc.cv.MB3, forks = ifelse(CV > 0, 'MB3_UF', 'MB3_LF'))

pdf("..", height = 4, width = 1.5)

ggplot(glyc.cv.MB3, aes(reorder(Gene, CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB3_UF", "MB3_LF")), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#762a83','#d7301f'))+
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  ylab(label = 'CV value') +
  theme(legend.position="none")

dev.off()

# TCA cycle ----

tca.df <- subset(cv.diff.metab, Gene %in% c('PDHA1','PDHA2','PDHB','PDHX','CS','ACO1','ACO2','IDH1','IDH2','IDH3A','IDH3B','IDH3G','OGDH','OGDHL','DLST','SUCLG1','SUCLG2','SDHA','SDHB','SDHC','SDHD','FH1','DHTKD1','L2HGDH','DLD','MDH1','MDH2','PCK1','PCK2','ACLY','DLAT','PC'))
tca.df <- tca.df[order(tca.df$MB1_CV),]
tca.df.g <- gather(tca.df, `MB1_CV`, `MB2_CV`, key = "bicluster", value = "CV")
tca.df.g.forks <- tca.df.g %>% mutate (forks = ifelse(bicluster == 'MB1_CV' & CV > '0', 'MB1_UF', 
                                                        ifelse(bicluster == 'MB1_CV' & CV < '0', 'MB1_LF', 
                                                               ifelse(bicluster == 'MB2_CV' & CV > '0', 'MB2_UF', 'MB2_LF'))))
tca.df.g.forks$MB1_CV <- as.numeric(paste(tca.df$MB1_CV)) 



pdf("..", height = 4, width = 3)


ggplot(tca.df.g.forks, aes(reorder(Gene, MB1_CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB2_UF", "MB1_UF", 'MB2_LF', 'MB1_LF')), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#1b7837','#D55E00','#980043','#a6bddb'))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  # scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2))+
  ylab(label = 'CV value') +
  theme(legend.position="none") +
  ylim(-1.3, 2)

dev.off()

# including MB3, and create a similar plot only with that bicluster:

cv.MB3 <- data.frame(df.list$MB3.all.metab.genes.cv)

tca.cv.MB3 <- subset(cv.MB3, Gene %in% c('PDHA1','PDHA2','PDHB','PDHX','CS','ACO1','ACO2','IDH1','IDH2','IDH3A','IDH3B','IDH3G','OGDH','OGDHL','DLST','SUCLG1','SUCLG2','SDHA','SDHB','SDHC','SDHD','FH1','DHTKD1','L2HGDH','DLD','MDH1','MDH2','PCK1','PCK2','ACLY','DLAT','PC'))
tca.cv.MB3 <- tca.cv.MB3[order(tca.cv.MB3$CV),]
tca.cv.MB3 <- mutate(tca.cv.MB3, forks = ifelse(CV > 0, 'MB3_UF', 'MB3_LF'))

pdf("..", height = 4, width = 1.5)

ggplot(tca.cv.MB3, aes(reorder(Gene, CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB3_UF", "MB3_LF")), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#762a83','#d7301f'))+
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  ylab(label = 'CV value') +
  theme(legend.position="none")

dev.off()


# pyruvate metabolism ----

pyr.df <- subset(cv.diff.metab, Gene %in% c('PCX','PCK2','MDH1','ME2','ME3','MOD1','ME1','ME3','PDHA1','PDHA2','PDHB','PDHX','DLAT','MDH1','MDH2','PC','PKLR','PKM2','PCK1','PCK2','LDHA','LDHAL6B','LDHB','LDHC','LDHD','ACSS1','ACSS2','ACACA','ACACB','ACAT1','ACOT12','ALDH1B1','ALDH2','ALDH3A2','ALDH7A1','ALDH9A1'))
pyr.df <- pyr.df[order(pyr.df$MB1_CV),]
pyr.df.g <- gather(pyr.df, `MB1_CV`, `MB2_CV`, key = "bicluster", value = "CV")
pyr.df.g.forks <- pyr.df.g %>% mutate (forks = ifelse(bicluster == 'MB1_CV' & CV > '0', 'MB1_UF', 
                                                      ifelse(bicluster == 'MB1_CV' & CV < '0', 'MB1_LF', 
                                                             ifelse(bicluster == 'MB2_CV' & CV > '0', 'MB2_UF', 'MB2_LF'))))
pyr.df.g.forks$MB1_CV <- as.numeric(paste(pyr.df$MB1_CV)) 



pdf("..", height = 4, width = 3)


ggplot(pyr.df.g.forks, aes(reorder(Gene, MB1_CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB2_UF", "MB1_UF", 'MB2_LF', 'MB1_LF')), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#1b7837','#D55E00','#980043','#a6bddb'))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  # scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2))+
  ylab(label = 'CV value') +
  theme(legend.position="none") +
  ylim(-1.3, 2)

dev.off()

# including MB3, and create a similar plot only with that bicluster:

cv.MB3 <- data.frame(df.list$MB3.all.metab.genes.cv)

pyr.cv.MB3 <- subset(cv.MB3, Gene %in% c('PCX','PCK2','MDH1','ME2','ME3','MOD1','ME1','ME3','PDHA1','PDHA2','PDHB','PDHX','DLAT','MDH1','MDH2','PC','PKLR','PKM2','PCK1','PCK2','LDHA','LDHAL6B','LDHB','LDHC','LDHD','ACSS1','ACSS2','ACACA','ACACB','ACAT1','ACOT12','ALDH1B1','ALDH2','ALDH3A2','ALDH7A1','ALDH9A1'))
pyr.cv.MB3 <- pyr.cv.MB3[order(pyr.cv.MB3$CV),]
pyr.cv.MB3 <- mutate(pyr.cv.MB3, forks = ifelse(CV > 0, 'MB3_UF', 'MB3_LF'))

pdf("..", height = 4, width = 1.5)

ggplot(pyr.cv.MB3, aes(reorder(Gene, CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB3_UF", "MB3_LF")), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#762a83','#d7301f'))+
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  ylab(label = 'CV value') +
  theme(legend.position="none")

dev.off()


# pentose phosphate pathway ----

ppp.df <- subset(cv.diff.metab, Gene %in% c('GPI', 'RBKS', 'PGM2', 'DERA', 'PRPS2','TKTL1','RPIA','RPE','ALDOA','PFKM','FBP1','PGD','PGLS','G6PD','IDNK','RGN','TALDO1','PGM1','GLYCTK','H6PD','PGD'))
ppp.df <- ppp.df[order(ppp.df$MB1_CV),]
ppp.df.g <- gather(ppp.df, `MB1_CV`, `MB2_CV`, key = "bicluster", value = "CV")
ppp.df.g.forks <- ppp.df.g %>% mutate (forks = ifelse(bicluster == 'MB1_CV' & CV > '0', 'MB1_UF', 
                                                      ifelse(bicluster == 'MB1_CV' & CV < '0', 'MB1_LF', 
                                                             ifelse(bicluster == 'MB2_CV' & CV > '0', 'MB2_UF', 'MB2_LF'))))
ppp.df.g.forks$MB1_CV <- as.numeric(paste(ppp.df$MB1_CV)) 



pdf("..", height = 4, width = 3)


ggplot(ppp.df.g.forks, aes(reorder(Gene, MB1_CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB2_UF", "MB1_UF", 'MB2_LF', 'MB1_LF')), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#1b7837','#D55E00','#980043','#a6bddb'))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  # scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2))+
  ylab(label = 'CV value') +
  theme(legend.position="none") +
  ylim(-1.3, 2)

dev.off()

# fatty acid biosynthesis ----

fab.df <- subset(cv.diff.metab, Gene %in% c('FASN', 'OLAH', 'OXSM','MCAT','ACACA','ACSL1','ACSF3','MECR','CBR4','HSD17B8'))
fab.df <- fab.df[order(fab.df$MB1_CV),]
fab.df.g <- gather(fab.df, `MB1_CV`, `MB2_CV`, key = "bicluster", value = "CV")
fab.df.g.forks <- fab.df.g %>% mutate (forks = ifelse(bicluster == 'MB1_CV' & CV > '0', 'MB1_UF', 
                                                      ifelse(bicluster == 'MB1_CV' & CV < '0', 'MB1_LF', 
                                                             ifelse(bicluster == 'MB2_CV' & CV > '0', 'MB2_UF', 'MB2_LF'))))
fab.df.g.forks$MB1_CV <- as.numeric(paste(fab.df$MB1_CV)) 



pdf("..", height = 4, width = 3)


ggplot(fab.df.g.forks, aes(reorder(Gene, MB1_CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB2_UF", "MB1_UF", 'MB2_LF', 'MB1_LF')), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#1b7837','#D55E00','#980043','#a6bddb'))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  # scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2))+
  ylab(label = 'CV value') +
  theme(legend.position="none") +
  ylim(-1.3, 2)

dev.off()


# including MB3, and create a similar plot only with that bicluster:

cv.MB3 <- data.frame(df.list$MB3.all.metab.genes.cv)

fab.cv.MB3 <- subset(cv.MB3, Gene %in% c('FASN', 'OLAH', 'OXSM','MCAT','ACACA','ACSL1','ACSF3','MECR','CBR4','HSD17B8'))
fab.cv.MB3 <- fab.cv.MB3[order(fab.cv.MB3$CV),]
fab.cv.MB3 <- mutate(fab.cv.MB3, forks = ifelse(CV > 0, 'MB3_UF', 'MB3_LF'))

pdf("..", height = 4, width = 1.5)

ggplot(fab.cv.MB3, aes(reorder(Gene, CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB3_UF", "MB3_LF")), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#762a83','#d7301f'))+
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  ylab(label = 'CV value') +
  theme(legend.position="none")

dev.off()


# fatty acid degradation ----

fad.df <- subset(cv.diff.metab, Gene %in% c('ACADS', 'HADH', 'ECHS1','ACADM','ECI1','ACAT2','HADHB', 'HADHA', 'ACADL', 'ADH1A','GCDH','CYP2U1','ACSL1','ACOX1','ACADVL','CPT1A','CPT2','ALDH1B1','ACADSB'))
fad.df <- fad.df[order(fad.df$MB1_CV),]
fad.df.g <- gather(fad.df, `MB1_CV`, `MB2_CV`, key = "bicluster", value = "CV")
fad.df.g.forks <- fad.df.g %>% mutate (forks = ifelse(bicluster == 'MB1_CV' & CV > '0', 'MB1_UF', 
                                                      ifelse(bicluster == 'MB1_CV' & CV < '0', 'MB1_LF', 
                                                             ifelse(bicluster == 'MB2_CV' & CV > '0', 'MB2_UF', 'MB2_LF'))))
fad.df.g.forks$MB1_CV <- as.numeric(paste(fad.df$MB1_CV)) 



pdf("..", height = 4, width = 3)


ggplot(fad.df.g.forks, aes(reorder(Gene, MB1_CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB2_UF", "MB1_UF", 'MB2_LF', 'MB1_LF')), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#1b7837','#D55E00','#980043','#a6bddb'))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  # scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2))+
  ylab(label = 'CV value') +
  theme(legend.position="none") +
  ylim(-1.3, 2)

dev.off()


cv.MB3 <- data.frame(df.list$MB3.all.metab.genes.cv)

fad.cv.MB3 <- subset(cv.MB3, Gene %in% c('ACADS', 'HADH', 'ECHS1','ACADM','ECI1','ACAT2','HADHB', 'HADHA', 'ACADL', 'ADH1A','GCDH','CYP2U1','ACSL1','ACOX1','ACADVL','CPT1A','CPT2','ALDH1B1','ACADSB'))
fad.cv.MB3 <- fad.cv.MB3[order(fad.cv.MB3$CV),]
fad.cv.MB3 <- mutate(fad.cv.MB3, forks = ifelse(CV > 0, 'MB3_UF', 'MB3_LF'))

pdf("..", height = 4, width = 1.5)

ggplot(fad.cv.MB3, aes(reorder(Gene, CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB3_UF", "MB3_LF")), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#762a83','#d7301f'))+
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  ylab(label = 'CV value') +
  theme(legend.position="none")

dev.off()




# purine metabolism ----

TEMP <- scan(text= 'RRM1 XDH PRPS2 PPAT GART ATIC ADSL PAICS PDE1A ADCY3 ENTPD3 ENTPD2 AK5 ADSS AMPD1 NT5E APRT ADA PNP NME1 GMPR IMPDH1 GMPS CANT1 HPRT1 GUK1 GDA ALLC HDDC3 PRUNE GUCY1A2 PAPSS2 PFAS ITPA NUDT2 DCK PGM2 NUDT9 FHIT NTPCR URAD PDE4A PGM1 NUDT16 AK3 ENPP3 HDDC2 ADK', what="")
dput(TEMP)

pur.df <- subset(cv.diff.metab, Gene %in% dput(TEMP))
pur.df <- pur.df[order(pur.df$MB1_CV),]
pur.df.g <- gather(pur.df, `MB1_CV`, `MB2_CV`, key = "bicluster", value = "CV")
pur.df.g.forks <- pur.df.g %>% mutate (forks = ifelse(bicluster == 'MB1_CV' & CV > '0', 'MB1_UF', 
                                                      ifelse(bicluster == 'MB1_CV' & CV < '0', 'MB1_LF', 
                                                             ifelse(bicluster == 'MB2_CV' & CV > '0', 'MB2_UF', 'MB2_LF'))))
pur.df.g.forks$MB1_CV <- as.numeric(paste(pur.df$MB1_CV)) 



pdf("..", height = 4, width = 3)


ggplot(pur.df.g.forks, aes(reorder(Gene, MB1_CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB2_UF", "MB1_UF", 'MB2_LF', 'MB1_LF')), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#1b7837','#D55E00','#980043','#a6bddb'))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  # scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2))+
  ylab(label = 'CV value') +
  theme(legend.position="none") +
  ylim(-1.3, 2)

dev.off()

# including MB3, and create a similar plot only with that bicluster:

cv.MB3 <- data.frame(df.list$MB3.all.metab.genes.cv)

pur.cv.MB3 <- subset(cv.MB3, Gene %in% c("RRM1", "XDH", "PRPS2", "PPAT", "GART", "ATIC", "ADSL", "PAICS", 
                                         "PDE1A", "ADCY3", "ENTPD3", "ENTPD2", "AK5", "ADSS", "AMPD1", 
                                         "NT5E", "APRT", "ADA", "PNP", "NME1", "GMPR", "IMPDH1", "GMPS", 
                                         "CANT1", "HPRT1", "GUK1", "GDA", "ALLC", "HDDC3", "PRUNE", "GUCY1A2", 
                                         "PAPSS2", "PFAS", "ITPA", "NUDT2", "DCK", "PGM2", "NUDT9", "FHIT", 
                                         "NTPCR", "URAD", "PDE4A", "PGM1", "NUDT16", "AK3", "ENPP3", "HDDC2", 
                                         "ADK"))
pur.cv.MB3 <- pur.cv.MB3[order(pur.cv.MB3$CV),]
pur.cv.MB3 <- mutate(pur.cv.MB3, forks = ifelse(CV > 0, 'MB3_UF', 'MB3_LF'))

pdf("..", height = 4, width = 1.5)

ggplot(pur.cv.MB3, aes(reorder(Gene, CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB3_UF", "MB3_LF")), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#762a83','#d7301f'))+
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  ylab(label = 'CV value') +
  theme(legend.position="none")

dev.off()


# pyrimidine metabolism ----

TEMP <- scan(text= 'RRM1
CAD
UMPS
NME1
CMPK1
CTPS2
ENTPD3
CANT1
UCK2
NT5E
UPP1
CDA
DCTPP1
DCK
DCTD
DTYMK
TK1
TYMP
TYMS
DPYD
DPYS
UPB1
DHODH
ENPP3
NUDT2
HDDC2
DUT
ASMTL
NT5C3A', what="")
dput(TEMP)

pyr.df <- subset(cv.diff.metab, Gene %in% dput(TEMP))
pyr.df <- pyr.df[order(pyr.df$MB1_CV),]
pyr.df.g <- gather(pyr.df, `MB1_CV`, `MB2_CV`, key = "bicluster", value = "CV")
pyr.df.g.forks <- pyr.df.g %>% mutate (forks = ifelse(bicluster == 'MB1_CV' & CV > '0', 'MB1_UF', 
                                                      ifelse(bicluster == 'MB1_CV' & CV < '0', 'MB1_LF', 
                                                             ifelse(bicluster == 'MB2_CV' & CV > '0', 'MB2_UF', 'MB2_LF'))))
pyr.df.g.forks$MB1_CV <- as.numeric(paste(pyr.df$MB1_CV)) 



pdf("..", height = 4, width = 3)


ggplot(pyr.df.g.forks, aes(reorder(Gene, MB1_CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB2_UF", "MB1_UF", 'MB2_LF', 'MB1_LF')), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#1b7837','#D55E00','#980043','#a6bddb'))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  # scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2))+
  ylab(label = 'CV value') +
  theme(legend.position="none") +
  ylim(-1.3, 2)

dev.off()

# including MB3, and create a similar plot only with that bicluster:

cv.MB3 <- data.frame(df.list$MB3.all.metab.genes.cv)

pyr.cv.MB3 <- subset(cv.MB3, Gene %in% c("RRM1", "CAD", "UMPS", "NME1", "CMPK1", "CTPS2", "ENTPD3", 
                                         "CANT1", "UCK2", "NT5E", "UPP1", "CDA", "DCTPP1", "DCK", "DCTD", 
                                         "DTYMK", "TK1", "TYMP", "TYMS", "DPYD", "DPYS", "UPB1", "DHODH", 
                                         "ENPP3", "NUDT2", "HDDC2", "DUT", "ASMTL", "NT5C3A"))
pyr.cv.MB3 <- pyr.cv.MB3[order(pyr.cv.MB3$CV),]
pyr.cv.MB3 <- mutate(pyr.cv.MB3, forks = ifelse(CV > 0, 'MB3_UF', 'MB3_LF'))

pdf("..", height = 4, width = 1.5)

ggplot(pyr.cv.MB3, aes(reorder(Gene, CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB3_UF", "MB3_LF")), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#762a83','#d7301f'))+
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  ylab(label = 'CV value') +
  theme(legend.position="none")

dev.off()

# oxidative phosphorylation ----

TEMP <- scan(text= 'COX4I2 NDUFC2 NDUFC1 NDUFB10 NDUFB9 NDUFB8 NDUFB7 NDUFB6 NDUFB5 NDUFB4 NDUFB3 NDUFB2 NDUFB1 NDUFA11 NDUFAB1 NDUFA10 NDUFA9 NDUFA8 NDUFA7 NDUFA6 NDUFA5 NDUFA4L2 NDUFA3 NDUFA2 NDUFA1 NDUFV3 NDUFV2 NDUFV1 NDUFS8 NDUFS7 NDUFS6 NDUFS5 NDUFS4 NDUFS3 NDUFS2 NDUFS1 ND6 ND5 ND4L ND4 ND3 ND2 ND1 COX17 COX15 COX11 COX8A COX7C COX7B COX7A1 COX6C CYC1 COX5B ATP5A1 ATP12A PPA1 SDHC SDHD SDHA SDHB UQCR11 UQCR10 UQCRQ UQCRB UQCRH UQCRC2 UQCRC1 CYTB UQCRFS1 COX6B1 COX6A1 COX5A COX2 COX10 COX1 COX3 ATP5B ATP5C1 ATP5O ATP5D ATP5E ATP5G2 ATP6 ATP5F1 ATP5I ATP5J ATP5J2 ATP8 ATP5H ATP5L ATP6V1A ATP6V1B1 ATP6V1C1 ATP6V1D ATP6V1E1 ATP6V1F ATP6V1G2 ATP6V1H ATP6V0A1 ATP6V0D1 ATP6V0E1 ATP6AP1 ATP6V0C NDUFB11 NDUFA12 NDUFA13 CYCS' , what="")

oxphos.df <- subset(cv.diff.metab, Gene %in% dput(TEMP))
oxphos.df <- oxphos.df[order(oxphos.df$MB1_CV),]
oxphos.df.g <- gather(oxphos.df, `MB1_CV`, `MB2_CV`, key = "bicluster", value = "CV")
oxphos.df.g.forks <- oxphos.df.g %>% mutate (forks = ifelse(bicluster == 'MB1_CV' & CV > '0', 'MB1_UF', 
                                                      ifelse(bicluster == 'MB1_CV' & CV < '0', 'MB1_LF', 
                                                             ifelse(bicluster == 'MB2_CV' & CV > '0', 'MB2_UF', 'MB2_LF'))))
oxphos.df.g.forks$MB1_CV <- as.numeric(paste(oxphos.df$MB1_CV)) 



pdf("..", height = 4, width = 3)


ggplot(oxphos.df.g.forks, aes(reorder(Gene, MB1_CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB2_UF", "MB1_UF", 'MB2_LF', 'MB1_LF')), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#1b7837','#D55E00','#980043','#a6bddb'))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 4),
        axis.title.x = element_text(size = 6, vjust = 2),
        axis.text.x = element_text(size = 6),
        axis.line.x = element_line(size = 0.2),
        axis.ticks.x = element_line(size = 0.2),
        axis.ticks.length.x = unit(0.1, 'cm')
  ) +
  # scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2))+
  ylab(label = 'CV value') +
  theme(legend.position="none") +
  ylim(-1, 2)

dev.off()

# including MB3, and create a similar plot only with that bicluster:

cv.MB3 <- data.frame(df.list$MB3.all.metab.genes.cv)

oxphos.cv.MB3 <- subset(cv.MB3, Gene %in% c("COX4I2", "NDUFC2", "NDUFC1", "NDUFB10", "NDUFB9", "NDUFB8", 
                                            "NDUFB7", "NDUFB6", "NDUFB5", "NDUFB4", "NDUFB3", "NDUFB2", "NDUFB1", 
                                            "NDUFA11", "NDUFAB1", "NDUFA10", "NDUFA9", "NDUFA8", "NDUFA7", 
                                            "NDUFA6", "NDUFA5", "NDUFA4L2", "NDUFA3", "NDUFA2", "NDUFA1", 
                                            "NDUFV3", "NDUFV2", "NDUFV1", "NDUFS8", "NDUFS7", "NDUFS6", "NDUFS5", 
                                            "NDUFS4", "NDUFS3", "NDUFS2", "NDUFS1", "ND6", "ND5", "ND4L", 
                                            "ND4", "ND3", "ND2", "ND1", "COX17", "COX15", "COX11", "COX8A", 
                                            "COX7C", "COX7B", "COX7A1", "COX6C", "CYC1", "COX5B", "ATP5A1", 
                                            "ATP12A", "PPA1", "SDHC", "SDHD", "SDHA", "SDHB", "UQCR11", "UQCR10", 
                                            "UQCRQ", "UQCRB", "UQCRH", "UQCRC2", "UQCRC1", "CYTB", "UQCRFS1", 
                                            "COX6B1", "COX6A1", "COX5A", "COX2", "COX10", "COX1", "COX3", 
                                            "ATP5B", "ATP5C1", "ATP5O", "ATP5D", "ATP5E", "ATP5G2", "ATP6", 
                                            "ATP5F1", "ATP5I", "ATP5J", "ATP5J2", "ATP8", "ATP5H", "ATP5L", 
                                            "ATP6V1A", "ATP6V1B1", "ATP6V1C1", "ATP6V1D", "ATP6V1E1", "ATP6V1F", 
                                            "ATP6V1G2", "ATP6V1H", "ATP6V0A1", "ATP6V0D1", "ATP6V0E1", "ATP6AP1", 
                                            "ATP6V0C", "NDUFB11", "NDUFA12", "NDUFA13", "CYCS"))
oxphos.cv.MB3 <- oxphos.cv.MB3[order(oxphos.cv.MB3$CV),]
oxphos.cv.MB3 <- mutate(oxphos.cv.MB3, forks = ifelse(CV > 0, 'MB3_UF', 'MB3_LF'))

pdf("..", height = 4, width = 1.5)

ggplot(oxphos.cv.MB3, aes(reorder(Gene, CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB3_UF", "MB3_LF")), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#762a83','#d7301f'))+
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 4),
        axis.title.x = element_text(size = 6, vjust = 2),
        axis.text.x = element_text(size = 6),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  ylab(label = 'CV value') +
  theme(legend.position="none")

dev.off()

# mitochondrial DNA maintenance ----

cv.diff.mito <- data.frame(df.list$MB1.allgenome.cv, df.list$MB2.allgenome.cv[,2], ((df.list$MB1.allgenome.cv[,2]) - (df.list$MB2.allgenome.cv[,2])), stringsAsFactors = FALSE)
colnames(cv.diff.mito) <- c('Gene', 'MB1_CV', 'MB2_CV', 'MB1-MB2')



TEMP <- scan(text= 'APEX1 ATAD3A ATAD3B DNA2 ENDOG EXOG LIG3 METTL4 MGME1 MTERF1 MTERF2 MUTYH OGG1 PIF1 POLB POLDIP2 POLG POLG2 POLQ POLRMT PPA2 PRIMPOL RECQL4 RNASEH1 SSBP1 TFAM TFB2M TOP1MT TOP3A TWNK UNG' , what="")

mtDNA.df <- subset(cv.diff.mito, Gene %in% dput(TEMP))
mtDNA.df <- mtDNA.df[order(mtDNA.df$MB1_CV),]
mtDNA.df.g <- gather(mtDNA.df, `MB1_CV`, `MB2_CV`, key = "bicluster", value = "CV")
mtDNA.df.g.forks <- mtDNA.df.g %>% mutate (forks = ifelse(bicluster == 'MB1_CV' & CV > '0', 'MB1_UF', 
                                                            ifelse(bicluster == 'MB1_CV' & CV < '0', 'MB1_LF', 
                                                                   ifelse(bicluster == 'MB2_CV' & CV > '0', 'MB2_UF', 'MB2_LF'))))
mtDNA.df.g.forks$MB1_CV <- as.numeric(paste(mtDNA.df$MB1_CV)) 



pdf("..", height = 4, width = 3)


ggplot(mtDNA.df.g.forks, aes(reorder(Gene, MB1_CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB2_UF", "MB1_UF", 'MB2_LF', 'MB1_LF')), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#1b7837','#D55E00','#980043','#a6bddb'))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  # scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2))+
  ylab(label = 'CV value') +
  theme(legend.position="none") +
  ylim(-1, 2)

dev.off()

# including MB3, and create a similar plot only with that bicluster:

cv.MB3 <- data.frame(df.list$MB3.allgenome.cv)

mtDNA.cv.MB3 <- subset(cv.MB3, Gene %in% c("APEX1", "ATAD3A", "ATAD3B", "DNA2", "ENDOG", "EXOG", "LIG3", 
                                           "METTL4", "MGME1", "MTERF1", "MTERF2", "MUTYH", "OGG1", "PIF1", 
                                           "POLB", "POLDIP2", "POLG", "POLG2", "POLQ", "POLRMT", "PPA2", 
                                           "PRIMPOL", "RECQL4", "RNASEH1", "SSBP1", "TFAM", "TFB2M", "TOP1MT", 
                                           "TOP3A", "TWNK", "UNG"))
mtDNA.cv.MB3 <- mtDNA.cv.MB3[order(mtDNA.cv.MB3$CV),]
mtDNA.cv.MB3 <- mutate(mtDNA.cv.MB3, forks = ifelse(CV > 0, 'MB3_UF', 'MB3_LF'))

pdf("..", height = 4, width = 1.5)

ggplot(mtDNA.cv.MB3, aes(reorder(Gene, CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB3_UF", "MB3_LF")), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#762a83','#d7301f'))+
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  ylab(label = 'CV value') +
  theme(legend.position="none")

dev.off()


# mitochondrial RNA metabolism ----

cv.diff.mito <- data.frame(df.list$MB1.allgenome.cv, df.list$MB2.allgenome.cv[,2], ((df.list$MB1.allgenome.cv[,2]) - (df.list$MB2.allgenome.cv[,2])), stringsAsFactors = FALSE)
colnames(cv.diff.mito) <- c('Gene', 'MB1_CV', 'MB2_CV', 'MB1-MB2')



TEMP <- scan(text= 'ALKBH1 ANGEL2 CDK5RAP1 DDX28 DHX30 DUS2 ELAC2 ENDOG ERAL1 FASTK FASTKD1 FASTKD2 FASTKD3 FASTKD5 GRSF1 GTPBP3 HSD17B10 LACTB2 LRPPRC METTL15 METTL17 METTL5 METTL8 MRM1 MRM2 MRM3 MRPL12 MRPL47 MRPS7 MRPS9 MTERF1 MTO1 MTPAP MTRES1 MYG1 NGRN NOA1 NSUN2 NSUN3 NSUN4 OSGEPL1 PDE12 PNPT1 POLRMT PPA2 PRORP PTCD1 PTCD2 PUS1 QTRT1 RCC1L REXO2 RMND1 RNASEH1 RPUSD3 RPUSD4 SLIRP SUPV3L1 TBRG4 TEFM TFAM TFB1M TFB2M THG1L TOP1MT TRIT1 TRMT1 TRMT10C TRMT2B TRMT5 TRMT61B TRMU TRNT1 TRUB2 YBEY YRDC' , what="")

mtRNA.df <- subset(cv.diff.mito, Gene %in% dput(TEMP))
mtRNA.df <- mtRNA.df[order(mtRNA.df$MB1_CV),]
mtRNA.df.g <- gather(mtRNA.df, `MB1_CV`, `MB2_CV`, key = "bicluster", value = "CV")
mtRNA.df.g.forks <- mtRNA.df.g %>% mutate (forks = ifelse(bicluster == 'MB1_CV' & CV > '0', 'MB1_UF', 
                                                          ifelse(bicluster == 'MB1_CV' & CV < '0', 'MB1_LF', 
                                                                 ifelse(bicluster == 'MB2_CV' & CV > '0', 'MB2_UF', 'MB2_LF'))))
mtRNA.df.g.forks$MB1_CV <- as.numeric(paste(mtRNA.df$MB1_CV)) 



pdf("..", height = 4, width = 3)


ggplot(mtRNA.df.g.forks, aes(reorder(Gene, MB1_CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB2_UF", "MB1_UF", 'MB2_LF', 'MB1_LF')), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#1b7837','#D55E00','#980043','#a6bddb'))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  # scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2))+
  ylab(label = 'CV value') +
  theme(legend.position="none") +
  ylim(-1.3, 2)

dev.off()

# including MB3, and create a similar plot only with that bicluster:

cv.MB3 <- data.frame(df.list$MB3.allgenome.cv)

mtRNA.cv.MB3 <- subset(cv.MB3, Gene %in% c("ALKBH1", "ANGEL2", "CDK5RAP1", "DDX28", "DHX30", "DUS2", "ELAC2", 
                                           "ENDOG", "ERAL1", "FASTK", "FASTKD1", "FASTKD2", "FASTKD3", "FASTKD5", 
                                           "GRSF1", "GTPBP3", "HSD17B10", "LACTB2", "LRPPRC", "METTL15", 
                                           "METTL17", "METTL5", "METTL8", "MRM1", "MRM2", "MRM3", "MRPL12", 
                                           "MRPL47", "MRPS7", "MRPS9", "MTERF1", "MTO1", "MTPAP", "MTRES1", 
                                           "MYG1", "NGRN", "NOA1", "NSUN2", "NSUN3", "NSUN4", "OSGEPL1", 
                                           "PDE12", "PNPT1", "POLRMT", "PPA2", "PRORP", "PTCD1", "PTCD2", 
                                           "PUS1", "QTRT1", "RCC1L", "REXO2", "RMND1", "RNASEH1", "RPUSD3", 
                                           "RPUSD4", "SLIRP", "SUPV3L1", "TBRG4", "TEFM", "TFAM", "TFB1M", 
                                           "TFB2M", "THG1L", "TOP1MT", "TRIT1", "TRMT1", "TRMT10C", "TRMT2B", 
                                           "TRMT5", "TRMT61B", "TRMU", "TRNT1", "TRUB2", "YBEY", "YRDC"))
mtRNA.cv.MB3 <- mtRNA.cv.MB3[order(mtRNA.cv.MB3$CV),]
mtRNA.cv.MB3 <- mutate(mtRNA.cv.MB3, forks = ifelse(CV > 0, 'MB3_UF', 'MB3_LF'))

pdf("..", height = 4, width = 1.5)

ggplot(mtRNA.cv.MB3, aes(reorder(Gene, CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB3_UF", "MB3_LF")), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#762a83','#d7301f'))+
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 9, vjust = 2),
        axis.text.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.3),
        axis.ticks.x = element_line(size = 0.3),
        axis.ticks.length.x = unit(0.15, 'cm')
  ) +
  ylab(label = 'CV value') +
  theme(legend.position="none")

dev.off()


# mitochondrial translation ----

cv.diff.mito <- data.frame(df.list$MB1.allgenome.cv, df.list$MB2.allgenome.cv[,2], ((df.list$MB1.allgenome.cv[,2]) - (df.list$MB2.allgenome.cv[,2])), stringsAsFactors = FALSE)
colnames(cv.diff.mito) <- c('Gene', 'MB1_CV', 'MB2_CV', 'MB1-MB2')



TEMP <- scan(text= 'AARS2 AURKAIP1 C12orf65 CARS2 CHCHD1 COA3 COX14 DAP3 DARS2 DDX28 DHX30 EARS2 ERAL1 EXD2 FARS2 FASTKD2 GADD45GIP1 GARS1 GATB GATC GFM1 GFM2 GRSF1 GTPBP10 GUF1 HARS2 HEMK1 IARS2 KARS1 LARS2 LRPPRC MALSU1 MARS2 METAP1D METTL17 MIEF1 MPV17L2 MRM2 MRM3 MRPL1 MRPL10 MRPL11 MRPL12 MRPL13 MRPL14 MRPL15 MRPL16 MRPL17 MRPL18 MRPL19 MRPL2 MRPL20 MRPL21 MRPL22 MRPL23 MRPL24 MRPL27 MRPL28 MRPL3 MRPL30 MRPL32 MRPL33 MRPL34 MRPL35 MRPL36 MRPL37 MRPL38 MRPL39 MRPL4 MRPL40 MRPL41 MRPL42 MRPL43 MRPL44 MRPL45 MRPL46 MRPL47 MRPL48 MRPL49 MRPL50 MRPL51 MRPL52 MRPL53 MRPL54 MRPL55 MRPL57 MRPL58 MRPL9 MRPS10 MRPS11 MRPS12 MRPS14 MRPS15 MRPS16 MRPS17 MRPS18A MRPS18B MRPS18C MRPS2 MRPS21 MRPS22 MRPS23 MRPS24 MRPS25 MRPS26 MRPS27 MRPS28 MRPS30 MRPS31 MRPS33 MRPS34 MRPS35 MRPS36 MRPS5 MRPS6 MRPS7 MRPS9 MRRF MTERF3 MTERF4 MTFMT MTG1 MTG2 MTIF2 MTIF3 MTRES1 MTRF1 MTRF1L NARS2 NGRN NOA1 NSUN4 OXA1L PARS2 PDF PPA2 PTCD3 PUSL1 QRSL1 RARS2 RBFA RMND1 SARS2 SLIRP TACO1 TARS2 TFB1M TIMM21 TRMT61B TSFM TUFM VARS2 WARS2 YARS2 YBEY' , what="")

mt_ribo.df <- subset(cv.diff.mito, Gene %in% dput(TEMP))
mt_ribo.df <- mt_ribo.df[order(mt_ribo.df$MB1_CV),]
mt_ribo.df.g <- gather(mt_ribo.df, `MB1_CV`, `MB2_CV`, key = "bicluster", value = "CV")
mt_ribo.df.g.forks <- mt_ribo.df.g %>% mutate (forks = ifelse(bicluster == 'MB1_CV' & CV > '0', 'MB1_UF', 
                                                          ifelse(bicluster == 'MB1_CV' & CV < '0', 'MB1_LF', 
                                                                 ifelse(bicluster == 'MB2_CV' & CV > '0', 'MB2_UF', 'MB2_LF'))))
mt_ribo.df.g.forks$MB1_CV <- as.numeric(paste(mt_ribo.df$MB1_CV)) 



pdf("..", height = 4, width = 3)


ggplot(mt_ribo.df.g.forks, aes(reorder(Gene, MB1_CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB2_UF", "MB1_UF", 'MB2_LF', 'MB1_LF')), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#1b7837','#D55E00','#980043','#a6bddb'))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 2),
        axis.title.x = element_text(size = 4, vjust = 2),
        axis.text.x = element_text(size = 4),
        axis.line.x = element_line(size = 0.15),
        axis.ticks.x = element_line(size = 0.15),
        axis.ticks.length.x = unit(0.1, 'cm')
  ) +
  # scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2))+
  ylab(label = 'CV value') +
  theme(legend.position="none") +
  ylim(-1, 2)

dev.off()

# including MB3, and create a similar plot only with that bicluster:

cv.MB3 <- data.frame(df.list$MB3.allgenome.cv)

mt_ribo.cv.MB3 <- subset(cv.MB3, Gene %in% c("AARS2", "AURKAIP1", "C12orf65", "CARS2", "CHCHD1", "COA3", 
                                             "COX14", "DAP3", "DARS2", "DDX28", "DHX30", "EARS2", "ERAL1", 
                                             "EXD2", "FARS2", "FASTKD2", "GADD45GIP1", "GARS1", "GATB", "GATC", 
                                             "GFM1", "GFM2", "GRSF1", "GTPBP10", "GUF1", "HARS2", "HEMK1", 
                                             "IARS2", "KARS1", "LARS2", "LRPPRC", "MALSU1", "MARS2", "METAP1D", 
                                             "METTL17", "MIEF1", "MPV17L2", "MRM2", "MRM3", "MRPL1", "MRPL10", 
                                             "MRPL11", "MRPL12", "MRPL13", "MRPL14", "MRPL15", "MRPL16", "MRPL17", 
                                             "MRPL18", "MRPL19", "MRPL2", "MRPL20", "MRPL21", "MRPL22", "MRPL23", 
                                             "MRPL24", "MRPL27", "MRPL28", "MRPL3", "MRPL30", "MRPL32", "MRPL33", 
                                             "MRPL34", "MRPL35", "MRPL36", "MRPL37", "MRPL38", "MRPL39", "MRPL4", 
                                             "MRPL40", "MRPL41", "MRPL42", "MRPL43", "MRPL44", "MRPL45", "MRPL46", 
                                             "MRPL47", "MRPL48", "MRPL49", "MRPL50", "MRPL51", "MRPL52", "MRPL53", 
                                             "MRPL54", "MRPL55", "MRPL57", "MRPL58", "MRPL9", "MRPS10", "MRPS11", 
                                             "MRPS12", "MRPS14", "MRPS15", "MRPS16", "MRPS17", "MRPS18A", 
                                             "MRPS18B", "MRPS18C", "MRPS2", "MRPS21", "MRPS22", "MRPS23", 
                                             "MRPS24", "MRPS25", "MRPS26", "MRPS27", "MRPS28", "MRPS30", "MRPS31", 
                                             "MRPS33", "MRPS34", "MRPS35", "MRPS36", "MRPS5", "MRPS6", "MRPS7", 
                                             "MRPS9", "MRRF", "MTERF3", "MTERF4", "MTFMT", "MTG1", "MTG2", 
                                             "MTIF2", "MTIF3", "MTRES1", "MTRF1", "MTRF1L", "NARS2", "NGRN", 
                                             "NOA1", "NSUN4", "OXA1L", "PARS2", "PDF", "PPA2", "PTCD3", "PUSL1", 
                                             "QRSL1", "RARS2", "RBFA", "RMND1", "SARS2", "SLIRP", "TACO1", 
                                             "TARS2", "TFB1M", "TIMM21", "TRMT61B", "TSFM", "TUFM", "VARS2", 
                                             "WARS2", "YARS2", "YBEY"))
mt_ribo.cv.MB3 <- mt_ribo.cv.MB3[order(mt_ribo.cv.MB3$CV),]
mt_ribo.cv.MB3 <- mutate(mt_ribo.cv.MB3, forks = ifelse(CV > 0, 'MB3_UF', 'MB3_LF'))

pdf("..", height = 4, width = 1.5)

ggplot(mt_ribo.cv.MB3, aes(reorder(Gene, CV, decreasing = FALSE), CV, fill = factor(forks, levels=c("MB3_UF", "MB3_LF")), position = 'stack')) +
  geom_col(aes(alpha=abs(CV))) +
  scale_fill_manual(values = c('#762a83','#d7301f'))+
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1))+
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 2),
        axis.title.x = element_text(size = 4, vjust = 2),
        axis.text.x = element_text(size = 4),
        axis.line.x = element_line(size = 0.15),
        axis.ticks.x = element_line(size = 0.15),
        axis.ticks.length.x = unit(0.1, 'cm')
  ) +
  ylab(label = 'CV value') +
  theme(legend.position="none")

dev.off()




