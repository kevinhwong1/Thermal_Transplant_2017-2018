#Title: Thermal Transplant 2017-2018 - Surface Area (Adult 2017)
#Author: KH Wong
#Date Last Modified: 20210118
#See Readme file for details 

library(dbplyr)
library(tidyverse)
library(readr)
library(stringr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(Rmisc)

# Import Data
A2017.SA <- read.csv("data/2017/Surface.Area/Colony.2017.sa.csv")


Box.A2017.SA <- ggplot(A2017.SA, aes(x=Origin, y=sa, fill = Treatment)) +
  geom_boxplot(width=.5, outlier.colour=NA, position = position_dodge(width = 0.9)) +
  geom_point(pch = 21, position=position_jitterdodge(dodge.width=0.9), outlier.shape= NA) +
  #  stat_summary(fun.y=median, geom="line", position = position_dodge(width = 0.9), aes(group=Parental.Treatment))  + 
  #  stat_summary(fun.y=median, geom="point", position = position_dodge(width = 0.9)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
  scale_fill_manual(values=c("#FFFFFF", "#999999")) +
  xlab("Origin") + ylab(expression("Surface Area " (cm^{-2}))) + #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black")) +
  theme(legend.position = c(0.8,0.92), 
        legend.title = element_text(colour="black", size=26),
        legend.text = element_text(colour="black", size = 24))

mean.col.size <- summarySE(A2017.SA, measurevar="sa", groupvars=c("Origin"))

mean.col.size.treat <- summarySE(A2017.SA, measurevar="sa", groupvars=c("Origin", "Treatment"))
