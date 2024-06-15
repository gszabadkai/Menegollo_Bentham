
load("../METABRIC_starting_data_final_corr_groups.Rdata")

library(dplyr)
library(ggplot2)
library(ggthemes)
library(tidyr)
library(scales)



all.clinical.df2.g <- gather(all.clinical.df2, 'MB1.forkscale.fork_0.2', 'MB2.forkscale.fork_0.2', 'MB3.forkscale.fork_0.2', key = 'forkscales', value = 'forks')



# plot theme----

theme_fork3 <- theme_tufte(base_size = 16, base_family = 'sans') +
  theme(aspect.ratio = 1/1) +
  theme(axis.ticks.length = unit(8, "pt")) +
  theme(axis.ticks = element_line(colour = "black", size = 0.6)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 2, r = 0, b = 0, l = 0))) +
  theme(legend.text = element_text(size = 12)) +
  theme(axis.line = element_blank())+
  theme(axis.text.y = element_text(color = "black", size = 16))+
  theme(axis.text.x = element_text(color = "black", size = 16))+
  theme(legend.margin=margin(0,0,0,0), legend.box.margin=margin(20,0,0,-275), legend.justification="top")+
  theme(panel.grid.minor = element_blank())+
  theme(legend.title=element_text(size=14))+
  theme(legend.title.align = 0.5) +
  theme(plot.title = element_text(size = 15, face = 'plain', hjust = 0.1, vjust = -20))

theme_set(theme_fork3)
grids <- data.frame(c(0,0), c(-10,10))
colnames(grids) <- c('x','y')


# MB1:MB2----
# Fig.1D
pdf("..", height = 5.08, width = 10)

ggplot(all.clinical.df2.g %>% arrange(desc(forks)) %>% subset(forks %in% c('MB1_LF','MB1_UF', 'MB2_UF','MB2_LF')), aes(MB1.forkscale.log, MB2.forkscale.log)) +
  geom_point(data = all.clinical.df2.g, aes(MB1.forkscale.log, MB2.forkscale.log), colour = '#d9d9d9', alpha=0.1, shape = 16) +
  geom_point(aes(colour = forks), alpha=0.3, size=2.5, shape = 16) +
  scale_colour_manual(name = "Forks",values = c('MB1_UF' = '#D55E00', 
                                 'MB1_LF' = '#a6bddb', 
                                 'wMB1_none' = '#bdbdbd',
                                 'wMB2_none' = '#bdbdbd', 
                                 'MB2_UF' = '#1b7837', 
                                 'MB2_LF' = '#980043' )) +
  geom_rangeframe(data=data.frame(MB1.forkscale.log=c(-10,10), MB2.forkscale.log=c(-10,10)), size = 1.2) +
  scale_x_continuous(breaks = seq(-10, 10, 5), limits = c(-10,10), minor_breaks = seq(-20,20,1)) +
  scale_y_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.8))) +
  geom_line(data = grids, aes(x=x, y=y), linetype = "dotted") +
  geom_line(data = grids, aes(x=y, y=x), linetype = "dotted") +
  labs(x ="MB1 forkscale", y = "MB2 forkscale")
  
dev.off()



# PAM50 status
# PAM50 separate plots (Fig 3B)
pdf("..", height = 5.08, width = 10)

ggplot(all.clinical.df2.g %>% arrange(desc(forks)) %>% subset(PAM50 %in% c("Basal")), aes(MB1.forkscale.log, MB2.forkscale.log)) +
  geom_point(data = all.clinical.df2.g, aes(MB1.forkscale.log, MB2.forkscale.log), colour = '#969696', size=2.5, alpha=0.1, shape = 16) +
  geom_point(aes(colour = PAM50), alpha=0.2, size=2.5, shape = 16) +
  scale_color_manual(name = 'PAM50', values = c("Basal" = "#e41a1c", "Her2" = "#377eb8", "LumA"="#fec44f", "LumB"="#4daf4a", "Normal"="#bdbdbd"))+
  geom_rangeframe(data=data.frame(MB1.forkscale.log=c(-10,10), MB2.forkscale.log=c(-10,10)), size = 1.2) +
  scale_x_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  scale_y_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.8))) +
  geom_line(data = grids, aes(x=x, y=y), linetype = "dotted") +
  geom_line(data = grids, aes(x=y, y=x), linetype = "dotted") +
  labs(x ="MB1 forkscale", y = "MB2 forkscale")

dev.off()

pdf("..", height = 5.08, width = 10)

ggplot(all.clinical.df2.g %>% arrange(desc(forks)) %>% subset(PAM50 %in% c("LumA", "Normal")), aes(MB1.forkscale.log, MB2.forkscale.log)) +
  geom_point(data = all.clinical.df2.g, aes(MB1.forkscale.log, MB2.forkscale.log), colour = '#969696', size=2.5, alpha=0.1, shape = 16) +
  geom_point(aes(colour = PAM50), alpha=0.2, size=2.5, shape = 16) +
  scale_color_manual(name = 'PAM50', values = c("Basal" = "#e41a1c", "Her2" = "#377eb8", "LumA"="#fec44f", "LumB"="#4daf4a", "Normal"="#252525"))+
  geom_rangeframe(data=data.frame(MB1.forkscale.log=c(-10,10), MB2.forkscale.log=c(-10,10)), size = 1.2) +
  scale_x_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  scale_y_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.8))) +
  geom_line(data = grids, aes(x=x, y=y), linetype = "dotted") +
  geom_line(data = grids, aes(x=y, y=x), linetype = "dotted") +
  labs(x ="MB1 forkscale", y = "MB2 forkscale")

dev.off()

pdf("..", height = 5.08, width = 10)

ggplot(all.clinical.df2.g %>% arrange(desc(forks)) %>% subset(PAM50 %in% c('LumB',"Her2")), aes(MB1.forkscale.log, MB2.forkscale.log)) +
  geom_point(data = all.clinical.df2.g, aes(MB1.forkscale.log, MB2.forkscale.log), colour = '#969696', size=2.5, alpha=0.1, shape = 16) +
  geom_point(aes(colour = PAM50), alpha=0.6, size=2.5, shape = 16) +
  scale_color_manual(name = 'PAM50', values = c("Basal" = "#e41a1c", "Her2" = "#9ecae1", "LumA"="#fec44f", "LumB"="#4daf4a", "Normal"="#252525"))+
  geom_rangeframe(data=data.frame(MB1.forkscale.log=c(-10,10), MB2.forkscale.log=c(-10,10)), size = 1.2) +
  scale_x_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  scale_y_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.8))) +
  geom_line(data = grids, aes(x=x, y=y), linetype = "dotted") +
  geom_line(data = grids, aes(x=y, y=x), linetype = "dotted") +
  labs(x ="MB1 forkscale", y = "MB2 forkscale")

dev.off()


# MaSC scores (Fig3B)

pdf("..", height = 5.08, width = 10)

ggplot(all.clinical.df2.g %>% arrange(desc(forks)) %>% subset(PAM50 %in% c("Basal", "Her2", "LumA", "LumB", "Normal")), aes(MB1.forkscale.log, MB2.forkscale.log)) +
  geom_point(aes(colour = MaSC_sum), alpha=0.8, size=2.5, shape = 16) +
  scale_colour_gradientn(colours=c("#253494","#d9d9d9","#bd0026"), 
                         values = rescale(c(-0.4, 0, 0.6)), 
                         breaks = c(-0.3, 0, 0.3, 0.6), 
                         guide=guide_colourbar(title = "MaSC scale", title.hjust=0.5, title.vjust=2, barwidth = 0.7, barheight = 4, nbin=100, ticks = TRUE, label.theme = element_text(size = 10))) +
  geom_rangeframe(data=data.frame(MB1.forkscale.log=c(-10,10), MB2.forkscale.log=c(-10,10)), size = 1.2) +
  scale_x_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  scale_y_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  # guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.8))) +
  geom_line(data = grids, aes(x=x, y=y), linetype = "dotted") +
  geom_line(data = grids, aes(x=y, y=x), linetype = "dotted") +
  labs(x ="MB1 forkscale", y = "MB2 forkscale")

dev.off()

# LP scores (Fig3B)

pdf("..", height = 5.08, width = 10)

ggplot(all.clinical.df2.g %>% arrange(desc(forks)) %>% subset(PAM50 %in% c("Basal", "Her2", "LumA", "LumB", "Normal")), aes(MB1.forkscale.log, MB2.forkscale.log)) +
  geom_point(aes(colour = LP_sum), alpha=0.8, size=2.5, shape = 16) +
  scale_colour_gradientn(colours=c("#253494","#d9d9d9","#bd0026"), 
                         values = rescale(c(-0.45, 0, 0.55)), 
                         breaks = c(-0.4, 0, 0.4), 
                         guide=guide_colourbar(title = "LP scale", title.hjust=0.5, title.vjust=2, barwidth = 0.7, barheight = 4, nbin=100, ticks = TRUE, label.theme = element_text(size = 10))) +
  geom_rangeframe(data=data.frame(MB1.forkscale.log=c(-10,10), MB2.forkscale.log=c(-10,10)), size = 1.2) +
  scale_x_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  scale_y_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  # guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.8))) +
  geom_line(data = grids, aes(x=x, y=y), linetype = "dotted") +
  geom_line(data = grids, aes(x=y, y=x), linetype = "dotted") +
  labs(x ="MB1 forkscale", y = "MB2 forkscale")

dev.off()


# mL scores (Fig3B)
pdf("..", height = 5.08, width = 10)


ggplot(all.clinical.df2.g %>% arrange(desc(forks)) %>% subset(PAM50 %in% c("Basal", "Her2", "LumA", "LumB", "Normal")), aes(MB1.forkscale.log, MB2.forkscale.log)) +
  geom_point(aes(colour = mL_sum), alpha=0.8, size=2.5, shape = 16) +
  scale_colour_gradientn(colours=c("#253494","#d9d9d9","#bd0026"), 
                         values = rescale(c(-0.5, 0, 0.45)), 
                         breaks = c(-0.4, 0, 0.4), 
                         guide=guide_colourbar(title = "mL scale", title.hjust=0.5, title.vjust=2, barwidth = 0.7, barheight = 4, nbin=100, ticks = TRUE, label.theme = element_text(size = 10))) +
  geom_rangeframe(data=data.frame(MB1.forkscale.log=c(-10,10), MB2.forkscale.log=c(-10,10)), size = 1.2) +
  scale_x_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  scale_y_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  # guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.8))) +
  geom_line(data = grids, aes(x=x, y=y), linetype = "dotted") +
  geom_line(data = grids, aes(x=y, y=x), linetype = "dotted") +
  labs(x ="MB1 forkscale", y = "MB2 forkscale")

dev.off()


# MB1:MB3---- Fig.1E

# Forks:
pdf("..", height = 5.08, width = 10)

ggplot(all.clinical.df2.g %>% arrange(desc(forks)) %>% subset(forks %in% c('MB1_LF','MB1_UF', 'MB3_UF','MB3_LF')), aes(MB1.forkscale.log, MB3.forkscale.log)) +
  geom_point(aes(colour = forks), alpha=3/5, size=2.5, shape = 16) +
  scale_colour_manual(name = "Forks",values = c('MB1_UF' = '#D55E00', 
                                                'MB1_LF' = '#a6bddb', 
                                                'wMB1_none' = '#bdbdbd',
                                                'wMB3_none' = '#bdbdbd', 
                                                'MB3_UF' = '#762a83', 
                                                'MB3_LF' = '#d7301f' )) +
  geom_rangeframe(data=data.frame(MB1.forkscale.log=c(-10,10), MB3.forkscale.log=c(-10,10)), size = 1.2) +
  scale_x_continuous(breaks = seq(-10, 10, 5), limits = c(-10,10), minor_breaks = seq(-20,20,1)) +
  scale_y_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.8))) +
  geom_line(data = grids, aes(x=x, y=y), linetype = "dotted") +
  geom_line(data = grids, aes(x=y, y=x), linetype = "dotted") +
  labs(x ="MB1 forkscale", y = "MB3 forkscale")

dev.off()


# MB2:MB3---- Fig.1E

pdf("..", height = 5.08, width = 10)

ggplot(all.clinical.df2.g %>% arrange(desc(forks)) %>% subset(forks %in% c('MB2_LF','MB2_UF', 'MB3_UF','MB3_LF')), aes(MB3.forkscale.log, MB2.forkscale.log)) +
  geom_point(aes(colour = forks), alpha=3/5, size=2.5, shape = 16) +
  scale_colour_manual(name = "Forks",values = c('MB2_UF' = '#1b7837', 
                                                'MB2_LF' = '#980043', 
                                                'wMB2_none' = '#bdbdbd',
                                                'wMB3_none' = '#bdbdbd', 
                                                'MB3_UF' = '#762a83', 
                                                'MB3_LF' = '#d7301f'  )) +
  geom_rangeframe(data=data.frame(MB2.forkscale.log=c(-10,10), MB3.forkscale.log=c(-10,10)), size = 1.2) +
  scale_x_continuous(breaks = seq(-10, 10, 5), limits = c(-10,10), minor_breaks = seq(-20,20,1)) +
  scale_y_continuous(breaks = seq(-10, 10, 5), limits = c(-10, 10), minor_breaks = seq(-20,20,1)) +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.8))) +
  geom_line(data = grids, aes(x=x, y=y), linetype = "dotted") +
  geom_line(data = grids, aes(x=y, y=x), linetype = "dotted") +
  labs(x ="MB3 forkscale", y = "MB2 forkscale")

dev.off()



