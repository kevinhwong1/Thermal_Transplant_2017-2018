#Title: Thermal Transplant 2017-2018 - Cholorphyll A (Adult and Larvae)
#Author: KH Wong
#Date Last Modified: 20190320
#See Readme file for details 

library(dbplyr)
library(tidyverse)
library(readr)
library(stringr)

### Adult 2017 Chl A ###

#Load data 

A2017.chla.run1 <- read.csv("data/2017/Chlorophyll/Chl_20190117_Run1.csv")
A2017.chla.run2 <- read.csv("data/2017/Chlorophyll/Chl_20190117_Run2.csv")
A2017.chla.run3 <- read.csv("data/2017/Chlorophyll/Chl_20190117_Run3.csv")
A2017.chla.run4 <- read.csv("data/2017/Chlorophyll/Chl_20190117_Run4.csv")
A2017.chla.run5 <- read.csv("data/2017/Chlorophyll/Chl_20190117_Run5.csv")

A2017.vial.meta <- read.csv("data/2017/Metadata/BIOS2017_Adult_Vial.csv") #Vial metadata
A2017.frag <- read.csv("data/2017/Surface.Area/Adult.Frag.2017.Calculated.csv") #Fragment metadata
A2017.meta <- read.csv("data/2017/Metadata/Colony.Metadata.csv") #Metadata
A2017.well.chla <- read.csv("data/2017/Chlorophyll/Adult_Ch1_Metadata.csv")

# Adding Run number to datasets

A2017.chla.run1$Run <- 1
A2017.chla.run2$Run <- 2
A2017.chla.run3$Run <- 3
A2017.chla.run4$Run <- 4
A2017.chla.run5$Run <- 5

# Combining Datasets
A2017.chla.all <- rbind(A2017.chla.run1, A2017.chla.run2, A2017.chla.run3, A2017.chla.run4, A2017.chla.run5)

# Make a unique column 
A2017.chla.all$run.well <- paste(A2017.chla.all$Run, A2017.chla.all$Well, sep = "-")
A2017.well.chla$run.well <- paste(A2017.well.chla$Run, A2017.well.chla$Well, sep = "-")

# Attaching vial metadata
A2017.chla.meta <- merge(A2017.chla.all, A2017.well.chla, by = "run.well" )

#Removing Blank values
A2017.chla.data <- A2017.chla.meta %>%
  filter(Vial != "BLANK")

# Attaching colony metadata
A2017.chla.data.vial <- merge(A2017.chla.data, A2017.vial.meta, by = "Vial")
A2017.chla.data.vial <- A2017.chla.data.vial %>%
  select(`Chl.630`, `Chl.663`, `Chl.750`, coral.id, timepoint)

A2017.chla.data.vial$coral.tp <- paste(A2017.chla.data.vial$coral.id, A2017.chla.data.vial$timepoint)

A2017.frag$coral.tp <- paste(A2017.frag$coral.id, A2017.frag$Timepoint)
A2017.frag <- A2017.frag %>%
  select(coral.tp, Surface.Area, homo.vol)

A2017.chla.data.meta <- merge(A2017.chla.data.vial, A2017.frag, by = "coral.tp")

#Subtracting 750 (blank) from 630 and 633 values
A2017.chla.data.meta$abs.630.corr <- A2017.chla.data.meta$`Chl.630` - A2017.chla.data.meta$`Chl.750`
A2017.chla.data.meta$abs.663.corr <- A2017.chla.data.meta$`Chl.663` - A2017.chla.data.meta$`Chl.750`

#Chlorphyll A concentration equation 
A2017.chla.data.meta$chlA.ug.sample <- 11.43*A2017.chla.data.meta$abs.663.corr - 0.64*A2017.chla.data.meta$abs.630.corr

#Standardization
A2017.chla.data.meta$ChlA.ugcm2 <- (A2017.chla.data.meta$chlA.ug.sample * (1000/500) * A2017.chla.data.meta$homo.vol)/A2017.chla.data.meta$Surface.Area #Calculating concentration

# Merging reefzone data
A2017.chla.final <- merge(A2017.chla.data.meta, A2017.meta, by ="coral.id")

# Summarizing 
mean.chla.A2017 <- summarySE(A2017.chla.final, measurevar="ChlA.ugcm2", groupvars=c("treatment", "reef.zone", "timepoint"))
mean.chla.A2017$reef.treatment <- paste(mean.chla.A2017$reef.zone, mean.chla.A2017$treatment)

#Reordering factors
mean.chla.A2017$timepoint <- factor(mean.chla.A2017$timepoint, levels = c("Pre-Treatment","Post-Treatment"))

#Plotting 
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Treatment"
Chla2017Adult <- ggplot(mean.chla.A2017, aes(x=timepoint, y=ChlA.ugcm2, group=reef.treatment)) + 
  geom_line(position=pd, color="black")+
  geom_errorbar(aes(ymin=ChlA.ugcm2-se, ymax=ChlA.ugcm2+se), width=.1, position=pd, color="black") + #Error bars
  ylim(1,6)+
  xlab("Timepoint") + ylab(expression("Chlorophyll a " (mu*g ~ cm^{-2}))) + #Axis titles
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

Chla2017Adult
ggsave(file = "output/Graphs/A2017.Chla.pdf", Chla2017Adult)

# Statistics

A2017.Chla.anova <- lm(ChlA.ugcm2~reef.zone*treatment*timepoint, data = A2017.chla.final)
qqnorm(resid(A2017.Chla.anova))
qqline(resid(A2017.Chla.anova))

boxplot(resid(A2017.Chla.anova)~A2017.chla.final$reef.zone)
boxplot(resid(A2017.Chla.anova)~A2017.chla.final$treatment)
boxplot(resid(A2017.Chla.anova)~A2017.chla.final$timepoint)
anova(A2017.Chla.anova)

capture.output(anova(A2017.Chla.anova), file = "output/Statistics/A2017.Chla.csv")
