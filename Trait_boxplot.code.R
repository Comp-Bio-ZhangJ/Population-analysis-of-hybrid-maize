library(data.table)
library("ape")
library(tidyverse)
library(ggplot2)
library(ggExtra)
library("grid")
library("gridExtra")
library(ggpubr)
library(reshape2)
library(ggprism)
library(patchwork)
library(magrittr)
library(rstatix)
library(RColorBrewer)
theme_set(theme_bw(16))
mycols1 <- colorRampPalette(brewer.pal(12,"Paired"))(12)[c(5,3,1)]

#1. trait data and plots.
raw.data <- read.table("C1M296443013.geno.trait.txt", header=T)
colnames(raw.data)[4] <- c("Genotype")
#Draw the plot.
plot.data <- melt(raw.data[,c(1,4,5:8)], id.vars = c("Strain", "Genotype"))
plot.data$Phe <- gsub("(\\w+)\\.\\w+", "\\1", plot.data[,3], perl=T)
plot.data$Type <- gsub("\\w+\\.(\\w+)", "\\1", plot.data[,3], perl=T)
plot.data <- plot.data[, c(1,2,5:6,4)]
colnames(plot.data)[5] <- c("Value")
head(plot.data, 3)
#   Strain Genotype Phe Type    Value
# 1 F01M01       TT  PH BLUE 273.5855
# 2 F01M02       TT  PH BLUE 270.9233
# 3 F01M03       CT  PH BLUE 273.3865

plot.data$Phe <- factor(plot.data$Phe, levels=c("PH", "EH"))
plot.data$Type <- factor(plot.data$Type, levels=c("BLUE", "MPH"))
plot.data$Genotype <- factor(plot.data$Genotype, levels=c("CC", "CT", "TT"))

#Pair-wise t-tests.
stat.frame <- plot.data %>%
  group_by(Phe, Type) %>%
  t_test(Value ~ Genotype, alternative="two.sided") %>%
  add_significance(p.col="p", cutpoints =c(0,0.01,0.05,1), symbols=c("**", "*", "")) %>%
  add_xy_position(x = "Genotype", fun="max", step.increase=0.015)
#Change to a pretty orders.
stat.frame$y.position <- stat.frame$y.position[c(2,3,1, 5,6,4, 8,9,7, 11,12,10)]
#Draw the figure.
fig.list <- list()
fig.list[[1]] <- ggplot(plot.data, aes(x=Genotype, y=Value)) +
  geom_boxplot(aes(fill=Genotype)) +
  scale_fill_manual(values=mycols1) +
  ggh4x::facet_grid2(Phe ~ Type, scales = "free_y", independent = "y") +
  labs(y = "Value", x="") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        axis.text.x=element_text(colour="black", size = 10, angle = 0, vjust = 0.5, hjust=0.5),
        axis.text.y=element_text(colour="black", size = 10),
        axis.title.x=element_text(colour="black", size = 10),
        axis.title.y=element_text(colour="black", size = 12), 
        legend.position="none") +
  stat_pvalue_manual(stat.frame, label = "{p.signif}")
pdf(c("C1M296443013.pretty.boxplot.pdf"),width=3.2, height=5)
grid.arrange(grobs=fig.list, ncol=1, nrow=1)
dev.off()
