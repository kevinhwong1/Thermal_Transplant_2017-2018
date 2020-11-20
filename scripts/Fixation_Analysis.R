#Title:Coral reproduction methods
#Author: KH Wong
#Date Last Modified: 20200713
#See Readme file for details 

### Larval Fixation ###

#Load Libaries
library(dbplyr)
library(tidyverse)
library(readr)
library(stringr)
library(Rmisc)

#Load Data
Prefix <- read.csv("data/2018/Larval.Size/Fixation_Analysis/PreFixationSizing.csv")
Postfix <- read.csv("data/2018/Larval.Size/Fixation_Analysis/PostFixationSizing.csv")

#Making a unique column
Prefix$coral.date <- paste(Prefix$coral.id, Prefix$Date)
Postfix$coral.date <- paste(Postfix$coral.id, Postfix$Date)


#Calculating means
prefix.mean <- summarySE(Prefix, measurevar="Volume", groupvars=c("coral.id","Date", "coral.date"))
postfix.mean <- summarySE(Postfix, measurevar="Volume", groupvars=c("coral.id","Date", "coral.date"))

#Adding Status column 
prefix.mean$Status <- as.factor("Live")
postfix.mean$Status <- as.factor("Fixed")

#Binding columns together
all.vol <- rbind(prefix.mean, postfix.mean)

#Boxplot
Box <- ggplot(all.vol, aes(x=Status, y=Volume, fill = Status)) +
  geom_boxplot(width=.3, outlier.colour=NA, position = position_dodge(width = 0.9)) +
  geom_jitter(position = position_jitter(width = 0.1), size = 4) +
  #  stat_summary(fun.y=median, geom="line", position = position_dodge(width = 0.9), aes(group=Parental.Treatment))  + 
  #  stat_summary(fun.y=median, geom="point", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("#FFFFFF", "#999999")) +
  xlab("Status") + ylab(expression("Larval Volume " (mm^{3}))) + #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black")) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/L2018.fixation.test.pdf", Box, width = 11, height = 11, units = c("in"))


# Statistics
Vol.anova <- lm(Volume~Status, data = all.vol)
qqnorm(resid(Vol.anova))
qqline(resid(Vol.anova))

boxplot(resid(Vol.anova)~all.vol$Status)

t.test(Volume~Status, data = all.vol)

capture.output(t.test(Volume~Status, data = all.vol), file = "output/Statistics/L2018.fixation.test.csv")


# 
# #Scatterplot
# fixation.plot<- ggplot(data = pre.post.size, aes(x=Volume.x, y=Volume.y))+
#   ylab("Pre-fixation Volume")+ xlab("Post-fixation Volume") + 
#   geom_point()+
#   geom_smooth(method = "lm") +
#   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# #creating a linear model
# lm.fix <- lm (Volume.x ~ Volume.y, data = pre.post.size) #creating a linear model off of standard curve
# lm.fix.summary <- summary(lm.fix) #creating a summary outfut for the linear model
# lm.fix.summary 
# 
# #calculating percent change 
# pre.post.size$p.change <- ((pre.post.size$Volume.y - pre.post.size$Volume.x)/pre.post.size$Volume.x) * 100 
# 
# #adding metadata
# pre.post.size$coral.id <- str_split_fixed(pre.post.size$coral.date, " ", 2) #creating colony id column
# colony.metadata <- read_csv("~/MyProjects/Reproduction_Methods/Data/Colony.Metadata.csv")
# 
# pre.post.size.data <- merge(pre.post.size, colony.metadata, by = "coral.id")
# 
# pchange.mean <- summarySE(pre.post.size.data, measurevar="p.change", groupvars="reef.zone")
# 
# pd <- position_dodge(0.1) # moves object .05 to the left and right
# pchange.plot<- ggplot(data = pchange.mean, aes(x=reef.zone, y=p.change))+
#   ylab("Percent change in Volume")+ xlab("Reef zone") + 
#   geom_point()+
#   geom_errorbar(aes(ymin=p.change-se, ymax=p.change+se), width=.1, position=pd, color="black") + #Error bars
#   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
# 
# 
