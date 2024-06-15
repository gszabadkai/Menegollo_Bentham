#figure 1B forks in ggplot



library(MCbiclust)
library(dplyr)
library(ggplot2)
library(ggthemes)

# Load data ----
load('../METABRIC_starting_data_final_corr_groups.Rdata')

# the theme ----
theme_fork2 <- theme_tufte(base_size = 16, base_family = 'sans') +
  theme(aspect.ratio = 2/3) +
  theme(axis.ticks.length = unit(8, "pt")) +
  theme(axis.ticks = element_line(colour = "black", size = 0.6)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.text = element_text(size = 12)) +
  theme(axis.line = element_blank())+
  theme(axis.text.y = element_text(color = "black", size = 16))+
  theme(axis.text.x = element_text(color = "black", size = 16))+
  theme(legend.margin=margin(0,0,0,0), legend.box.margin=margin(40,0,0,-20), legend.justification="top")+
  theme(panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size = 28, face = 'italic', hjust = 0.5, vjust = 2))

theme_set(theme_fork2)




# PLOTTING the MB1 BICLUSTER ----
# PAM50

forkbaseFig1B_MB1_PAM50 <- ggplot(subset(all.clinical.df, PAM50 %in% c("Basal", "Her2", "LumA", "LumB", "Normal")), aes(MB1.index, MB1.pc1)) +
  geom_point(aes(color=`PAM50`), alpha=3/5, size=3, shape = 16)+
  scale_color_manual(values = c("Basal" = "#e41a1c", "Her2" = "#377eb8", "LumA"="#fec44f", "LumB"="#4daf4a", "Normal"="#bdbdbd"))+
  geom_rangeframe(data=data.frame(MB1.pc1=c(-30,30), MB1.index=c(1,2000)), size = 1.2) +
  scale_x_continuous(breaks = seq(0, 2000, 500), limits = c(0, 2000), minor_breaks = seq(150,1050,300)) +
  scale_y_continuous(breaks = seq(-60, 60, 15), limits = c(-30, 30), minor_breaks = waiver()) +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.9)))+
  geom_line(aes(x=MB1.index, y=(log(MB1.index)/0.02)), linetype = "dotted")+
  geom_line(aes(x=MB1.index, y=(MB1.index* -0.02)), linetype = "dotted")+
  labs(title="MB1",x ="MB1 ranking", y = "MB1 pc1")

forkbaseFig1B_MB1_PAM50
ggsave("..pdf", width = 7.62, height = 5.08)


# PLOTTING the MB3 BICLUSTER ----
# PAM50
forkbaseFig1B_MB3_PAM50 <- ggplot(subset(all.clinical.df, PAM50 %in% c("Basal", "Her2", "LumA", "LumB", "Normal")), aes(MB3.index, MB3.pc1.rev)) +
  geom_point(aes(color=`PAM50`), alpha=3/5, size=3, shape=16)+
  scale_color_manual(values = c("Basal" = "#e41a1c", "Her2" = "#377eb8", "LumA"="#fec44f", "LumB"="#4daf4a", "Normal"="#bdbdbd"))+
  geom_rangeframe(data=data.frame(MB3.pc1.rev=c(-20,40), MB3.index=c(1,2000)), size = 1.2) +
  scale_x_continuous(breaks = seq(0, 2000, 500), limits = c(0, 2000), minor_breaks = seq(150,1050,300)) +
  scale_y_continuous(breaks = seq(-60, 60, 10), limits = c(-20, 40), minor_breaks = waiver()) +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.9)))+
  geom_line(aes(x=MB3.index, y=(MB3.index*0.02)), linetype = "dotted")+
  geom_line(aes(x=MB3.index, y=(MB3.index* -0.02)), linetype = "dotted") +
  labs(title="MB3",x ="MB3 ranking", y = "MB3 pc1")

forkbaseFig1A_MB3_PAM50
ggsave("..pdf", width = 7.62, height = 5.08)



# PLOTTING MB2 BICLUSTER ----
# PAM50
forkbaseFig1B_MB2_PAM50 <- ggplot(subset(all.clinical.df, PAM50 %in% c("Basal", "Her2", "LumA", "LumB", "Normal")), aes(MB2.index, MB2.pc1)) +
  geom_point(aes(color=`PAM50`), alpha=3/5, size=3, shape = 16)+
  scale_color_manual(values = c("Basal" = "#e41a1c", "Her2" = "#377eb8", "LumA"="#fec44f", "LumB"="#4daf4a", "Normal"="#bdbdbd"))+
  geom_rangeframe(data=data.frame(MB2.pc1=c(-30,30), MB2.index=c(1,2000)), size = 1.2) +
  scale_x_continuous(breaks = seq(0, 2000, 500), limits = c(0, 2000), minor_breaks = seq(150,1050,300)) +
  scale_y_continuous(breaks = seq(-60, 60, 15), limits = c(-30, 30), minor_breaks = waiver()) +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.9)))+
  geom_line(aes(x=MB2.index, y=(MB2.index*0.02)), linetype = "dotted")+
  geom_line(aes(x=MB2.index, y=(MB2.index* -0.02)), linetype = "dotted")+
  labs(title="MB2",x ="MB2 ranking", y = "MB2 pc1")

forkbaseFig1B_MB2_PAM50
ggsave("..pdf", width = 7.62, height = 5.08)

