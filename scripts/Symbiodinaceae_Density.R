#Title: Thermal Transplant 2017-2018 - Symbiodinaceae Densitiy (Adult and Larvae)
#Author: KH Wong
#Date Last Modified: 20190320
#See Readme file for details 

library(dbplyr)
library(tidyverse)
library(readr)
library(stringr)

### 2017 Adult Symbiodinaceae Densities ###

#Import Data
A2017.adult.Zoox <- read.csv("data/2017/Zoox/Adult_2017_Zoox.csv")

A2017.vial.meta <- read.csv("data/2017/Metadata/BIOS2017_Adult_Vial.csv") #Vial metadata
A2017.frag <- read.csv("data/2017/Surface.Area/Adult.Frag.2017.Calculated.csv") #Fragment metadata
A2017.meta <- read.csv("data/2017/Metadata/Colony.Metadata.csv") #Metadata

#Calculating cells/Larvae from cells/count
A2017.adult.Zoox$Average <- rowMeans(A2017.adult.Zoox[ ,7:12]) #averaging counts
A2017.adult.Zoox$Cells.mL <- (A2017.adult.Zoox$Average/A2017.adult.Zoox$squares.count)/0.0001 #number of cells per mL - volume of haemocytometer

#Removing metadata columns
A2017.adult.Zoox2 <- A2017.adult.Zoox %>%
  select("coral.id", "timepoint", "Cells.mL")

#Making a unique column to merge dataframes
A2017.frag$coral.tp <- paste(A2017.frag$coral.id, A2017.frag$Timepoint)
A2017.adult.Zoox2$coral.tp <- paste(A2017.adult.Zoox2$coral.id, A2017.adult.Zoox2$timepoint)

#Merging dataframes
A2017.zoox.colony <- merge(A2017.adult.Zoox2, A2017.meta, by = "coral.id")
A2017.zoox.meta <- merge(A2017.zoox.colony, A2017.frag, by = "coral.tp")

#Removing metadata columns
A2017.zoox.meta <- A2017.zoox.meta %>%
  select("coral.id.x","treatment", "reef.zone", "timepoint", "Cells.mL", "coral.tp", "homo.vol", "Surface.Area")

#Accounting for homogenate volume 
A2017.zoox.meta$Cells.mL.vol <- A2017.zoox.meta$Cells.mL*A2017.zoox.meta$homo.vol 

#Normalizing to fragment surface area
A2017.zoox.meta$Cells.cm2 <- A2017.zoox.meta$Cells.mL.vol/A2017.zoox.meta$Surface.Area 
A2017.zoox.meta$Cells.cm2.x6 <- A2017.zoox.meta$Cells.cm2 / 1000000

#Summarizing
A2017.zoox.mean <- summarySE(A2017.zoox.meta, measurevar="Cells.cm2.x6", groupvars=c("treatment", "reef.zone", "timepoint"))
A2017.zoox.mean$reef.treatment <- paste(A2017.zoox.mean$reef.zone, A2017.zoox.mean$treatment)

#Reordering factors
A2017.zoox.mean$timepoint <- factor(A2017.zoox.mean$timepoint, levels = c("Pre-Treatment","Post-Treatment"))

#Plotting 
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Treatment"
zoox2017Adult <- ggplot(A2017.zoox.mean, aes(x=timepoint, y=Cells.cm2.x6, group=reef.treatment)) + 
  geom_line(position=pd, color="black")+
  geom_errorbar(aes(ymin=Cells.cm2.x6-se, ymax=Cells.cm2.x6 +se), width=.1, position=pd, color="black") + #Error bars
  #  ylim(1,6)+
  xlab("Timepoint") + ylab(expression("Number of cells " (x10^{6} ~ cm^{-2})))+ #Axis titles
  geom_point(aes(shape=reef.zone, color=treatment), size=4, position=pd)+
  scale_shape_manual(values=c(16,17),
                     name = "Reef Zone")+
  scale_color_manual(values=c("dodgerblue3", "tomato1"),
                     name = "Treatment")+ #colour modification
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.title.x = element_blank()) +
  theme(legend.position = c(0.9,0.8))

zoox2017Adult

ggsave(file = "output/Graphs/A2017.Zoox.pdf", zoox2017Adult)

## Statistics
zoox.A2017.anova <- lm(Cells.cm2.x6~reef.zone*treatment*timepoint, data = A2017.zoox.meta)
qqnorm(resid(zoox.A2017.anova))
qqline(resid(zoox.A2017.anova))

boxplot(resid(zoox.A2017.anova)~A2017.zoox.meta$reef.zone)
boxplot(resid(zoox.A2017.anova)~A2017.zoox.meta$treatment)
boxplot(resid(zoox.A2017.anova)~A2017.zoox.meta$timepoint)
anova(zoox.A2017.anova)

capture.output(anova(zoox.A2017.anova), file = "output/Statistics/A2017.Zoox.csv")


### 2017 Larval Symbiodinaceae Densities ###

#Import Data
L2017.zoox <- read.csv("data/2017/Zoox/Larval_CellCounts.csv")

#Cleaning up data frame
L2017.zoox <- L2017.zoox[,colSums(is.na(L2017.zoox))<nrow(L2017.zoox)] #Removing columns with all NAs
L2017.zoox.clean <-L2017.zoox[complete.cases(L2017.zoox), ]

#Calculating cells/Larvae from cells/count
L2017.zoox.clean$Average <- rowMeans(L2017.zoox.clean[ ,12:17]) #averaging counts
L2017.zoox.clean$Cells.mL <- (L2017.zoox.clean$Average/L2017.zoox.clean$squares.count)/0.0001 #number of cells per mL - volume of haemocytometer
L2017.zoox.clean$Cells.mL.Diluted <- L2017.zoox.clean$Cells.mL/L2017.zoox.clean$Dilution.Factor #factoring in dilution
L2017.zoox.clean$Cells.larvae <- L2017.zoox.clean$Cells.mL.Diluted/L2017.zoox.clean$Sample.Size #diving by number of larvae per vial

#Remove metadata columns
L2017.zoox.clean2 <- L2017.zoox.clean %>% select(c(Vial, Larval.Release.Date, Colony, Cells.larvae))

#Renaming columns
colnames(L2017.zoox.clean)[colnames(L2017.zoox.clean)=='Colony']<-'coral.id' #changing column name 
colnames(L2017.zoox.clean)[colnames(L2017.zoox.clean)=='Larval.Release.Date']<-'date' #changing column name 

#Attaching metadata
L2017.zoox.meta <- merge(L2017.zoox.clean, A2017.meta, by = "coral.id")

#Subsetting patch
L2017.zoox.patch <- subset(L2017.zoox.meta, reef.zone=="Patch")
L2017.zoox.patch$Cells.larvae.x3 <-  L2017.zoox.patch$Cells.larvae / 1000

#Summarizing
L2017.zoox.patch.mean <- summarySE(L2017.zoox.patch, measurevar="Cells.larvae.x3", groupvars="treatment")

#Plotting reaction norm
pd <- position_dodge(0.1) # moves object .05 to the left and right
legend.title <- "Treatment"
#legend.title <- "Parental Treatment"
zoox.larv.2017<- ggplot(L2017.zoox.patch.mean, aes(x=treatment, y=Cells.larvae.x3)) + 
  geom_line(position=pd, color="black")+
  geom_errorbar(aes(ymin=Cells.larvae.x3-se, ymax=Cells.larvae.x3+se), width=.1, position=pd, color="black") + #Error bars
  ylim(5,8)+
  xlab("Treatment") + ylab(expression("Number of cells " (x10^{3} ~ Larva^{-1})))+ #Axis titles
  geom_point(aes(fill=treatment), color ="black", pch=21, size=4)+
  scale_fill_manual(legend.title, values=c("dodgerblue3", "tomato1"))+ #colour modification
  theme_bw() + theme(panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.title.x = element_blank()) 

ggsave(file = "output/Graphs/L2017.Zoox.Patch.pdf", zoox.larv.2017)

#Statistics
L2017.zoox.Patch.anova <- lm(Cells.larvae.x3~treatment, data = L2017.zoox.patch)
qqnorm(resid(L2017.zoox.Patch.anova))
qqline(resid(L2017.zoox.Patch.anova))

boxplot(resid(L2017.zoox.Patch.anova)~L2017.zoox.patch$treatment)

t.test(Cells.larvae.x3~treatment, data = L2017.zoox.patch)

capture.output(t.test(Cells.larvae.x3~treatment, data = L2017.zoox.patch), file = "output/Statistics/L2017.Zoox.Patch.csv")
