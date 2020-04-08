#Title: Thermal Transplant 2017-2018 - Cholorphyll A (Adult and Larvae)
#Author: KH Wong
#Date Last Modified: 20190527
#See Readme file for details 

library(dbplyr)
library(tidyverse)
library(readr)
library(stringr)
library(ggplot2)
library(Rmisc)
library(gdata)


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

#T2 analysis
A2017.chla.T2 <- A2017.chla.final %>%
  filter(timepoint == "Post-Treatment")

# Summarizing 
mean.A2017.col.chla.T2 <- summarySE(A2017.chla.T2, measurevar="ChlA.ugcm2", groupvars=c("coral.id", "treatment", "reef.zone"))
mean.A2017.chla.T2 <- summarySE(A2017.chla.T2, measurevar="ChlA.ugcm2", groupvars=c("treatment", "reef.zone"))

# Plotting
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Treatment"
Chla2017Adult.T2 <- ggplot(mean.A2017.chla.T2, aes(x=treatment, y=ChlA.ugcm2, group=reef.zone)) + 
  geom_line(position=pd, color="black", aes(linetype=reef.zone), size = 2)+
  geom_errorbar(aes(ymin=ChlA.ugcm2-se, ymax=ChlA.ugcm2+se), width=.1, size = 1, position=pd, color="black") + #Error bars
  ylim(0,5.0)+
  xlab("Timepoint") + ylab(expression("Chlorophyll a " (mu*g ~ cm^{-2})))+ #Axis titles
  geom_point(aes(fill=treatment, shape=reef.zone), size=14, position=pd, color = "black")+
  scale_shape_manual(values=c(21,24),
                     name = "Reef Zone")+
  scale_fill_manual(values=c("#FFFFFF", "#999999"),
                    name = "Treatment")+ #colour modification
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/A2017.Chla.T2.pdf", Chla2017Adult.T2, width = 11, height = 11, units = c("in"))

## Statistics
A2017.Chla.anova.T2 <- lm(ChlA.ugcm2~reef.zone*treatment, data = mean.A2017.col.chla.T2)
qqnorm(resid(A2017.Chla.anova.T2))
qqline(resid(A2017.Chla.anova.T2))

boxplot(resid(A2017.Chla.anova.T2)~mean.A2017.col.chla.T2$reef.zone)
boxplot(resid(A2017.Chla.anova.T2)~mean.A2017.col.chla.T2$treatment)

anova(A2017.Chla.anova.T2)

capture.output(anova(A2017.Chla.anova.T2), file = "output/Statistics/A2017.Chla.T2.csv")

#Violin plot
Violin.A2017.Chla <- ggplot(mean.A2017.col.chla.T2, aes(x=reef.zone, y=ChlA.ugcm2, fill = treatment)) +
  geom_violin(position = position_dodge(width = 0.9)) +
  geom_boxplot(width=.3, outlier.colour=NA, position = position_dodge(width = 0.9)) +
  #  stat_summary(fun.y=median, geom="line", position = position_dodge(width = 0.9), aes(group=Parental.Treatment))  + 
  #  stat_summary(fun.y=median, geom="point", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("#FFFFFF", "#999999")) +
  xlab("Reef Zone") + ylab(expression("Chlorophyll a " (ug ~ cm^{-2}))) + #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black")) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/A2017.Chla.violin.pdf", Violin.A2017.Chla, width = 11, height = 11, units = c("in"))




### Larval 2017 Chl a ###

#Load data 

L2017.chla.run1 <- read.csv("data/2017/Chlorophyll/Chl_Larvae_20190307_Run1.csv")
L2017.chla.run2 <- read.csv("data/2017/Chlorophyll/Chl_Larvae_20190307_Run2.csv")
L2017.chla.run3 <- read.csv("data/2017/Chlorophyll/Chl_Larvae_20190307_Run3.csv")
L2017.chla.run4 <- read.csv("data/2017/Chlorophyll/Chl_Larvae_20190307_Run4.csv")
L2017.chla.run5 <- read.csv("data/2017/Chlorophyll/Chl_Larvae_20190307_Run5.csv")
L2017.chla.run6 <- read.csv("data/2017/Chlorophyll/Chl_Larvae_20190307_Run6.csv")
L2017.chla.run7 <- read.csv("data/2017/Chlorophyll/Chl_Larvae_20190307_Run7.csv")
L2017.chla.run8 <- read.csv("data/2017/Chlorophyll/Chl_Larvae_20190307_Run8.csv")
L2017.chla.run9 <- read.csv("data/2017/Chlorophyll/Chl_Larvae_20190307_Run9.csv")
L2017.chla.run10 <- read.csv("data/2017/Chlorophyll/Chl_Larvae_20190307_Run10.csv")
L2017.chla.run11 <- read.csv("data/2017/Chlorophyll/Chl_Larvae_20190307_Run11.csv")

L2017.vial.meta <- read.csv("data/2017/Metadata/Vial.Metadata.csv") #Vial metadata
A2017.meta <- read.csv("data/2017/Metadata/Colony.Metadata.csv") #Metadata
L2017.well.chla <- read.csv("data/2017/Chlorophyll/Larvae_Ch1_Metadata.csv")

# Adding Run number to datasets

L2017.chla.run1$Run <- 1
L2017.chla.run2$Run <- 2
L2017.chla.run3$Run <- 3
L2017.chla.run4$Run <- 4
L2017.chla.run5$Run <- 5
L2017.chla.run6$Run <- 6
L2017.chla.run7$Run <- 7
L2017.chla.run8$Run <- 8
L2017.chla.run9$Run <- 9
L2017.chla.run10$Run <- 10
L2017.chla.run11$Run <- 11

# Combining Datasets
L2017.chla.all <- rbind(L2017.chla.run1, L2017.chla.run2, L2017.chla.run3, L2017.chla.run4, L2017.chla.run5, L2017.chla.run6, L2017.chla.run7, L2017.chla.run8, L2017.chla.run9, L2017.chla.run10, L2017.chla.run11)

# Make a unique column 
L2017.chla.all$run.well <- paste(L2017.chla.all$Run, L2017.chla.all$Well, sep = "-")
L2017.well.chla$run.well <- paste(L2017.well.chla$Run, L2017.well.chla$Well, sep = "-")

# Attaching vial metadata
L2017.chla.meta <- merge(L2017.chla.all, L2017.well.chla, by = "run.well" )

#Removing NA values
L2017.chla.data <- L2017.chla.meta %>%
  filter(Vial != "NA")
L2017.chla.data2 <- L2017.chla.data %>%
  filter(Chl.663 != "NA",Chl.630 != "NA",Chl.750 != "NA")

#Attaching colony metadata
L2017.chla.data.vial <- merge(L2017.chla.data2, L2017.vial.meta, by = "Vial")
L2017.chla.data.vial <- L2017.chla.data.vial %>%
  select(`Chl.630`, `Chl.663`, `Chl.750`, coral.id, Timepoint, Larval.Release.Date)

L2017.chla.data.vial$coral.tp <- paste(L2017.chla.data.vial$coral.id, L2017.chla.data.vial$timepoint)
L2017.chla.data.vial$'Chl.630' <- as.numeric(L2017.chla.data.vial$'Chl.630') #set Chl630 as numeric values

#Subtracting 750 (blank) from 630 and 633 values
L2017.chla.data.vial$abs.630.corr <- L2017.chla.data.vial$`Chl.630` - L2017.chla.data.vial$`Chl.750`
L2017.chla.data.vial$abs.663.corr <- L2017.chla.data.vial$`Chl.663` - L2017.chla.data.vial$`Chl.750`

#Chlorphyll A concentration equation 
L2017.chla.data.vial$chlA.ug.sample <- 11.43*L2017.chla.data.vial$abs.663.corr - 0.64*L2017.chla.data.vial$abs.630.corr

#Standardization
L2017.chla.data.vial$chlA.nglarv <- (L2017.chla.data.vial$chlA.ug.sample/20)*1000 #Calculating concentration. Dividing by 20 larvae per vial and converting ug to ng

# Merging reefzone data
L2017.chla.final <- merge(L2017.chla.data.vial, A2017.meta, by ="coral.id")

# Summarizing 
mean.chla.col.L2017 <- summarySE(L2017.chla.final, measurevar="chlA.nglarv", groupvars=c("treatment", "reef.zone", "coral.id"))

mean.chla.L2017 <- summarySE(L2017.chla.final, measurevar="chlA.nglarv", groupvars=c("treatment", "reef.zone", "Timepoint"))

#Reordering factors
mean.chla.L2017$Timepoint <- factor(mean.chla.L2017$Timepoint, levels = c("Pre-Treatment","Post-Treatment"))

mean.chla.L2017.patch <- mean.chla.L2017 %>% 
  filter(reef.zone == "Patch")

#Plotting 
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Treatment"
Chla2017Larvae <- ggplot(mean.chla.L2017.patch, aes(x=treatment, y=chlA.nglarv, group = reef.zone)) + 
  geom_line(position=pd, aes(linetype=reef.zone, color = treatment), size = 2)+
  geom_errorbar(aes(ymin=chlA.nglarv-se, ymax=chlA.nglarv+se), width=.1, position=pd, color="black") + #Error bars
  #ylim(1,6)+
  xlab("Treatment") + ylab(expression("Chlorophyll a " (ng~larva^{-1}))) + #Axis titles
  geom_point(aes(fill=treatment, shape=reef.zone), size=14, position=pd, color = "black")+
  scale_shape_manual(values=c(21,24),
                     name = "Reef Zone")+
  scale_fill_manual(values=c("#FFFFFF", "#999999"),
                    name = "Treatment")+ #colour modification
  scale_color_manual(values=c("black", "#999999"),
                     name = "Treatment")+
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")
  
Chla2017Larvae
ggsave(file = "output/Graphs/L2017.Chla.pdf", Chla2017Larvae, width = 11, height = 11, units = c("in"))

#Statistics
L2017.Chla.Patch.anova <- lm(chlA.nglarv~treatment, data = mean.chla.col.L2017)
qqnorm(resid(L2017.Chla.Patch.anova))
qqline(resid(L2017.Chla.Patch.anova))

boxplot(resid(L2017.Chla.Patch.anova)~mean.chla.col.L2017$treatment)

t.test(chlA.nglarv~treatment, data = mean.chla.col.L2017)

capture.output(t.test(chlA.nglarv~treatment, data = mean.chla.col.L2017), file = "output/Statistics/L2017.Chla.Patch.csv")


### Chla standardized by larval volume 

L2017.chla.final$Date.coral.ID <- paste(L2017.chla.final$Larval.Release.Date, L2017.chla.final$coral.id, sep = "-")

L.Size.2017.Patch <- read.csv("data/2017/Larval.Size/Patch_mean_Larvalsize.csv")

L2017.chla.final.size <- merge(L.Size.2017.Patch, L2017.chla.final, by = "Date.coral.ID")

L2017.chla.final.size$chla.ng.mm3 <- L2017.chla.final.size$chlA.nglarv/L2017.chla.final.size$Volume_1

#L2017.TP.Patch.meta.1 <- L2017.TP.Patch.meta[-c(43), ] #outlier removal

Violin.L2017.Chla <- ggplot(L2017.chla.final.size, aes(x=Treatment.1, y=chla.ng.mm3, fill = Treatment.1)) +
  geom_violin(position = position_dodge(width = 0.9)) +
  geom_boxplot(width=.3, outlier.colour=NA, position = position_dodge(width = 0.9)) +
  #  stat_summary(fun.y=median, geom="line", position = position_dodge(width = 0.9), aes(group=Parental.Treatment))  + 
  #  stat_summary(fun.y=median, geom="point", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("#FFFFFF", "#999999")) +
  xlab("Parental Treatment") + ylab(expression("Chlorophyll a " (ng ~ mm^{-3}))) + #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black")) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/L2017.Chla.Patch.size.pdf", Violin.L2017.Chla, width = 11, height = 11, units = c("in"))

capture.output(t.test(chla.ng.mm3~Treatment.1, data = L2017.chla.final.size), file = "output/Statistics/L2017.Chla.size.csv")



### Adult 2018 chl A ###

# Loading Data

A2018.chla.run1.1 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190501_Run1.csv")
A2018.chla.run2.1 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190501_Run2.csv")
A2018.chla.run3.1 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190501_Run3.csv")
A2018.chla.run4.1 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190501_Run4.csv")
A2018.chla.run5.1 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190501_Run5.csv")
A2018.chla.run6.1 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190501_Run6.csv")
A2018.chla.run7.1 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190501_Run7.csv")
A2018.chla.run8.1 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190501_Run8.csv")
A2018.chla.run9.1 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190501_Run9.csv")
A2018.chla.run10.1 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190501_Run10.csv")
A2018.chla.run11.1 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190501_Run11.csv")
A2018.chla.run12.1 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190501_Run12.csv")

A2018.chla.run1.2 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190502_Run1.csv")
A2018.chla.run2.2 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190502_Run2.csv")
A2018.chla.run3.2 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190502_Run3.csv")
A2018.chla.run4.2 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190502_Run4.csv")
A2018.chla.run5.2 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190502_Run5.csv")
A2018.chla.run6.2 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190502_Run6.csv")
A2018.chla.run7.2 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190502_Run7.csv")
A2018.chla.run8.2 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190502_Run8.csv")


A2018.vial.meta <- read.csv("data/2018/Metadata/BIOS2018_Adult_Vial.csv") #Vial metadata
A2018.frag <- read.csv("data/2018/Surface.Area/Adult.Frag.2018.Calculated.csv") #Fragment metadata
A2018.meta <- read.csv("data/2018/Metadata/Sample_Info_Transp.csv") #Metadata
A2018.well.chla.1 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190501_Wells.csv")
A2018.well.chla.2 <- read.csv("data/2018/Chlorophyll/Chl_Adults_20190502_Wells.csv")
A2018.vol <- read.csv("data/2018/Metadata/BIOS2018_Frag_Vol.csv")

#Adding run number to datasets
A2018.chla.run1.1$Run <- 1
A2018.chla.run2.1$Run <- 2
A2018.chla.run3.1$Run <- 3
A2018.chla.run4.1$Run <- 4
A2018.chla.run5.1$Run <- 5
A2018.chla.run6.1$Run <- 6
A2018.chla.run7.1$Run <- 7
A2018.chla.run8.1$Run <- 8
A2018.chla.run9.1$Run <- 9
A2018.chla.run10.1$Run <- 10
A2018.chla.run11.1$Run <- 11
A2018.chla.run12.1$Run <- 12

A2018.chla.run1.2$Run <- 1
A2018.chla.run2.2$Run <- 2
A2018.chla.run3.2$Run <- 3
A2018.chla.run4.2$Run <- 4
A2018.chla.run5.2$Run <- 5
A2018.chla.run6.2$Run <- 6
A2018.chla.run7.2$Run <- 7
A2018.chla.run8.2$Run <- 8

# Combining Datasets
A2018.chla.1 <- rbind(A2018.chla.run1.1,
                      A2018.chla.run2.1,
                      A2018.chla.run3.1,
                      A2018.chla.run4.1,
                      A2018.chla.run5.1,
                      A2018.chla.run6.1,
                      A2018.chla.run7.1,
                      A2018.chla.run8.1,
                      A2018.chla.run9.1,
                      A2018.chla.run10.1,
                      A2018.chla.run11.1,
                      A2018.chla.run12.1)

A2018.chla.2 <- rbind(A2018.chla.run1.2,
                      A2018.chla.run2.2,
                      A2018.chla.run3.2,
                      A2018.chla.run4.2,
                      A2018.chla.run5.2,
                      A2018.chla.run6.2,
                      A2018.chla.run7.2,
                      A2018.chla.run8.2)

# Make a unique column 
A2018.chla.1$run.well <- paste(A2018.chla.1$Run, A2018.chla.1$Well, sep = "-")
A2018.chla.2$run.well <- paste(A2018.chla.2$Run, A2018.chla.2$Well, sep = "-")


A2018.well.chla.1$run.well <- paste(A2018.well.chla.1$Run, A2018.well.chla.1$Well, sep = "-")
A2018.well.chla.2$run.well <- paste(A2018.well.chla.2$Run, A2018.well.chla.2$Well, sep = "-")

# Attaching vial metadata
A2018.chla.meta.1 <- merge(A2018.chla.1, A2018.well.chla.1, by = "run.well" )
A2018.chla.meta.2 <- merge(A2018.chla.2, A2018.well.chla.2, by = "run.well" )

# Merging datasets
A2018.chla.meta <- rbind(A2018.chla.meta.1,A2018.chla.meta.2)

#Removing NA values
A2018.chla.data <- A2018.chla.meta %>%
  filter(Vial != "NA")

# Attaching vial metadata
A2018.chla.data.vial <- merge(A2018.chla.data, A2018.vial.meta, by = "Vial")
A2018.chla.data.vial <- A2018.chla.data.vial %>%
  select(coral.id,`Chl.630`, `Chl.663`, `Chl.750`)

# Attaching SA data
colnames(A2018.frag)[which(names(A2018.frag) == "Colony")] <- "coral.id"
A2018.frag <- A2018.frag %>%
  select(coral.id, Surface.Area)
A2018.chla.data.sa <- merge(A2018.chla.data.vial, A2018.frag, by = "coral.id")

# Attaching homogenate volumes
A2018.vol <- A2018.vol %>%
  select(coral.id, Homo.Vol.mL)
A2018.chla.data.meta <- merge(A2018.chla.data.sa, A2018.vol, by = "coral.id")

#Subtracting 750 (blank) from 630 and 633 values
A2018.chla.data.meta$abs.630.corr <- A2018.chla.data.meta$`Chl.630` - A2018.chla.data.meta$`Chl.750`
A2018.chla.data.meta$abs.663.corr <- A2018.chla.data.meta$`Chl.663` - A2018.chla.data.meta$`Chl.750`

#Chlorphyll A concentration equation 
A2018.chla.data.meta$chlA.ug.sample <- 11.43*A2018.chla.data.meta$abs.663.corr - 0.64*A2018.chla.data.meta$abs.630.corr

#Standardization
A2018.chla.data.meta$ChlA.ugcm2 <- (A2018.chla.data.meta$chlA.ug.sample * (1000/500) * A2018.chla.data.meta$Homo.Vol.mL)/A2018.chla.data.meta$Surface.Area #Calculating concentration

# Merging reefzone data
A2018.chla.final <- merge(A2018.chla.data.meta, A2018.meta, by ="coral.id") %>%
  select(coral.id, ChlA.ugcm2, Origin, Treatment, Transplant.Site, Group)

# Summarizing 
mean.chla.col.A2018 <- summarySE(A2018.chla.final, measurevar="ChlA.ugcm2", groupvars=c("coral.id", "Origin", "Treatment", "Transplant.Site"))
mean.chla.col.A2018$orig.coral.id <- mean.chla.col.A2018$coral.id %>% str_replace("-A","") %>% str_replace("-B","")


mean.chla.A2018 <- summarySE(A2018.chla.final, measurevar="ChlA.ugcm2", groupvars=c("Origin", "Treatment", "Transplant.Site"))
mean.chla.A2018$reef.treatment <- paste(mean.chla.A2018$Origin, mean.chla.A2018$Treatment)

#Plotting 
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Treatment"
Chla2018Adult <- ggplot(mean.chla.A2018, aes(x=Transplant.Site, y=ChlA.ugcm2, group=reef.treatment)) + 
  geom_line(position=pd, aes(linetype=Origin, color = Treatment), size = 2)+
  geom_errorbar(aes(ymin=ChlA.ugcm2-se, ymax=ChlA.ugcm2+se), width=.1, position=pd, color="black") + #Error bars
  ylim(1.5,3)+
  xlab("Transplant Site") + ylab(expression("Chlorophyll a " (mu*g ~ cm^{-2}))) + #Axis titles
  geom_point(aes(fill=Treatment, shape=Origin), size=14, position=pd, color = "black")+
  scale_shape_manual(values=c(21,24),
                     name = "Reef Zone")+
  scale_fill_manual(values=c("#FFFFFF", "#999999"),
                    name = "Treatment")+ #colour modification
  scale_color_manual(values=c("black", "#999999"),
                     name = "Treatment")+
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

Chla2018Adult
ggsave(file = "output/Graphs/A2018.Chla.pdf", Chla2018Adult, width = 11, height = 11, units = c("in"))



## Statistics
Chla2018Adult.2018.anova <- lm(ChlA.ugcm2~Origin*Treatment*Transplant.Site, data = mean.chla.col.A2018)
qqnorm(resid(Chla2018Adult.2018.anova))
qqline(resid(Chla2018Adult.2018.anova)) 

boxplot(resid(Chla2018Adult.2018.anova)~mean.chla.col.A2018$Origin)
boxplot(resid(Chla2018Adult.2018.anova)~mean.chla.col.A2018$Treatment) #not normal
boxplot(resid(Chla2018Adult.2018.anova)~mean.chla.col.A2018$Transplant.Site)

anova(Chla2018Adult.2018.anova)

capture.output(anova(Chla2018Adult.2018.anova), file = "output/Statistics/A2018.Chla.csv")

#Nested model accounting for coral id
fit <- aov(ChlA.ugcm2~Origin*Treatment*Transplant.Site + Error(orig.coral.id), data=mean.chla.col.A2018)
summary(fit)

capture.output(summary(fit), file = "output/Statistics/A2018.Chla.nested.csv")

# TukeyHSD(fit, "Origin")
# 
# library(TukeyC)
# tuk <- TukeyHSD(mean.chla.col.A2018,
#              model = 'ChlA.ugcm2~Origin*Treatment*Transplant.Site + Error(orig.coral.id)',
#              error = 'orig.coral.id',
#              which = 'Transplant.Site',
#              fl1=1,
#              sig.level = 0.05)
# 
# summary(tuk)

### Larval Chl A 2018 ###

L2018.chla.run1.1 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190509_Run1.csv")
L2018.chla.run2.1 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190509_Run2.csv")
L2018.chla.run3.1 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190509_Run3.csv")
L2018.chla.run4.1 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190509_Run4.csv")
L2018.chla.run5.1 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190509_Run5.csv")
L2018.chla.run6.1 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190509_Run6.csv")
L2018.chla.run7.1 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190509_Run7.csv")
L2018.chla.run8.1 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190509_Run8.csv")
L2018.chla.run9.1 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190509_Run9.csv")
L2018.chla.run10.1 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190509_Run10.csv")
L2018.chla.run11.1 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190509_Run11.csv")
L2018.chla.run12.1 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190509_Run12.csv")

L2018.chla.run1.2 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190510_Run1.csv")
L2018.chla.run2.2 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190510_Run2.csv")
L2018.chla.run3.2 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190510_Run3.csv")
L2018.chla.run4.2 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190510_Run4.csv")
L2018.chla.run5.2 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190510_Run5.csv")
L2018.chla.run6.2 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190510_Run6.csv")
L2018.chla.run7.2 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190510_Run7.csv")
L2018.chla.run8.2 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190510_Run8.csv")
L2018.chla.run9.2 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190510_Run9.csv")


L2018.vial.meta <- read.csv("data/2018/Metadata/Larval_2018_vial_metadata.csv") #Vial metadata
A2018.meta <- read.csv("data/2018/Metadata/Sample_Info_Transp.csv") #Metadata
L2018.well.chla.1 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190509.csv")
L2018.well.chla.2 <- read.csv("data/2018/Chlorophyll/Chl_Larvae_20190510.csv")

#Adding run number to datasets
L2018.chla.run1.1$Run <- 1
L2018.chla.run2.1$Run <- 2
L2018.chla.run3.1$Run <- 3
L2018.chla.run4.1$Run <- 4
L2018.chla.run5.1$Run <- 5
L2018.chla.run6.1$Run <- 6
L2018.chla.run7.1$Run <- 7
L2018.chla.run8.1$Run <- 8
L2018.chla.run9.1$Run <- 9
L2018.chla.run10.1$Run <- 10
L2018.chla.run11.1$Run <- 11
L2018.chla.run12.1$Run <- 12

L2018.chla.run1.2$Run <- 1
L2018.chla.run2.2$Run <- 2
L2018.chla.run3.2$Run <- 3
L2018.chla.run4.2$Run <- 4
L2018.chla.run5.2$Run <- 5
L2018.chla.run6.2$Run <- 6
L2018.chla.run7.2$Run <- 7
L2018.chla.run8.2$Run <- 8
L2018.chla.run9.2$Run <- 9

# Combining Datasets
L2018.chla.1 <- rbind(L2018.chla.run1.1,
                      L2018.chla.run2.1,
                      L2018.chla.run3.1,
                      L2018.chla.run4.1,
                      L2018.chla.run5.1,
                      L2018.chla.run6.1,
                      L2018.chla.run7.1,
                      L2018.chla.run8.1,
                      L2018.chla.run9.1,
                      L2018.chla.run10.1,
                      L2018.chla.run11.1,
                      L2018.chla.run12.1)

L2018.chla.2 <- rbind(L2018.chla.run1.2, 
                      L2018.chla.run2.2,
                      L2018.chla.run3.2,
                      L2018.chla.run4.2,
                      L2018.chla.run5.2,
                      L2018.chla.run6.2,
                      L2018.chla.run7.2,
                      L2018.chla.run8.2,
                      L2018.chla.run9.2)

# Make a unique column 
L2018.chla.1$run.well <- paste(L2018.chla.1$Run, L2018.chla.1$Well, sep = "-")
L2018.chla.2$run.well <- paste(L2018.chla.2$Run, L2018.chla.2$Well, sep = "-")

L2018.well.chla.1$run.well <- paste(L2018.well.chla.1$Run, L2018.well.chla.1$Well, sep = "-")
L2018.well.chla.2$run.well <- paste(L2018.well.chla.2$Run, L2018.well.chla.2$Well, sep = "-")

# Attaching vial metadata
L2018.chla.meta.1 <- merge(L2018.chla.1, L2018.well.chla.1, by = "run.well" )
L2018.chla.meta.2 <- merge(L2018.chla.2, L2018.well.chla.2, by = "run.well" )

# Merging datasets
L2018.chla.meta <- rbind(L2018.chla.meta.1,L2018.chla.meta.2)

#Removing NA values
L2018.chla.data <- L2018.chla.meta %>%
  filter(Vial != "NA")

#Attaching vial metadata
L2018.chla.data.vial <- merge(L2018.chla.data, L2018.vial.meta, by = "Vial")
colnames(L2018.chla.data.vial)[colnames(L2018.chla.data.vial)=="Colony"] <- "coral.id"
L2018.chla.data.vial2 <- merge(L2018.chla.data.vial, A2018.meta, by = "coral.id")
L2018.chla.data.vial2 <- L2018.chla.data.vial2 %>%
  select(coral.id,Vial,Date.x,`Chl.630`, `Chl.663`, `Chl.750`, Origin, Treatment.y, Transplant.Site)

#Subtracting 750 (blank) from 630 and 633 values
L2018.chla.data.vial2$abs.630.corr <- L2018.chla.data.vial2$`Chl.630` - L2018.chla.data.vial2$`Chl.750`
L2018.chla.data.vial2$abs.663.corr <- L2018.chla.data.vial2$`Chl.663` - L2018.chla.data.vial2$`Chl.750`

#Chlorphyll A concentration equation 
L2018.chla.data.vial2$chlA.ug.sample <- 11.43*L2018.chla.data.vial2$abs.663.corr - 0.64*L2018.chla.data.vial2$abs.630.corr

#Standardization
L2018.chla.data.vial2$chlA.nglarv <- (L2018.chla.data.vial2$chlA.ug.sample/20)*1000 #Calculating concentration. Dividing by 20 larvae per vial and converting ug to ng
colnames(L2018.chla.data.vial2)[colnames(L2018.chla.data.vial2)=="Treatment.y"] <- "Treatment"

# Summarizing 

mean.chla.col.L2018 <- summarySE(L2018.chla.data.vial2, measurevar="chlA.nglarv", groupvars=c("coral.id","Origin", "Treatment", "Transplant.Site"))
mean.chla.col.L2018$orig.coral.id <- mean.chla.col.L2018$coral.id %>% str_replace("-A","") %>% str_replace("-B","")

mean.chla.L2018 <- summarySE(L2018.chla.data.vial2, measurevar="chlA.nglarv", groupvars=c("Origin", "Treatment", "Transplant.Site"))
mean.chla.L2018$reef.treatment <- paste(mean.chla.L2018$Origin, mean.chla.L2018$Treatment)

#Plotting 
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Treatment"
Chla2018Larvae <- ggplot(mean.chla.L2018, aes(x=Transplant.Site, y=chlA.nglarv, group=reef.treatment)) + 
  geom_line(position=pd, aes(linetype=Origin, color = Treatment), size = 2)+
  geom_errorbar(aes(ymin=chlA.nglarv-se, ymax=chlA.nglarv+se), width=.1, position=pd, color="black") + #Error bars
  ylim(8,30)+
  xlab("Transplant Site") + ylab(expression("Chlorophyll a " (ng  ~ Larva^{-1}))) + #Axis titles
  geom_point(aes(fill=Treatment, shape=Origin), size=14, position=pd, color = "black")+
  scale_shape_manual(values=c(21,24),
                     name = "Reef Zone")+
  scale_fill_manual(values=c("#FFFFFF", "#999999"),
                    name = "Treatment")+ #colour modification
  scale_color_manual(values=c("black", "#999999"),
                     name = "Treatment")+
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

Chla2018Larvae
ggsave(file = "output/Graphs/L2018.Chla.pdf", Chla2018Larvae, width = 11, height = 11, units = c("in"))

## Statistics
Chla2018larvae.2018.anova <- lm(chlA.nglarv~Origin*Treatment*Transplant.Site, data = mean.chla.col.L2018)
qqnorm(resid(Chla2018larvae.2018.anova))
qqline(resid(Chla2018larvae.2018.anova)) 

boxplot(resid(Chla2018larvae.2018.anova)~mean.chla.col.L2018$Origin)
boxplot(resid(Chla2018larvae.2018.anova)~mean.chla.col.L2018$Treatment) 
boxplot(resid(Chla2018larvae.2018.anova)~mean.chla.col.L2018$Transplant.Site)

anova(Chla2018larvae.2018.anova)

capture.output(anova(Chla2018larvae.2018.anova), file = "output/Statistics/L2018.Chla.csv")

# #Nested ANOVA
# Nested.L2018 <- aov(chlA.nglarv~Origin*Treatment*Transplant.Site + Error(orig.coral.id), data=mean.chla.col.L2018)
# summary(Nested.L2018) 

#### ng chla per mm3 ####

L2018.chla.data.vial2$Date.coral.ID <- paste(L2018.chla.data.vial2$Date.x, L2018.chla.data.vial2$coral.id, sep = "-")
#L2018.chla.colday <- summarySE(L2018.chla.data.vial2, measurevar="chlA.nglarv", groupvars=c("Date.coral.ID","Origin", "Treatment", "Transplant.Site"))

L.Size.2018 <- read.csv("data/2018/Larval.Size/Mean_Larvalsize_ColDay.csv")

L2018.Chla.meta.size <- merge(L.Size.2018, L2018.chla.data.vial2, by = "Date.coral.ID")

L2018.Chla.meta.size$Chla.ng.mm3 <- L2018.Chla.meta.size$chlA.nglarv/L2018.Chla.meta.size$Volume

# Summarizing 
L2018.Chla.mean.size <- summarySE(L2018.Chla.meta.size, measurevar="Chla.ng.mm3", groupvars=c("Origin", "Treatment", "Transplant.Site"))

#mean.TP.L2018 <- summarySE(mean.TP.col.L2018, measurevar="Conc.calcS", groupvars=c("Origin", "Treatment.y", "Transplant.Site"))
L2018.Chla.mean.size$reef.treatment <- paste(L2018.Chla.mean.size$Origin, L2018.Chla.mean.size$Treatment)

#Plotting 
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Treatment"
Chla2018Larvae.size <- ggplot(L2018.Chla.mean.size, aes(x=Transplant.Site, y=Chla.ng.mm3, group=reef.treatment)) + 
  geom_line(position=pd, aes(linetype=Origin, color = Treatment), size = 2)+
  geom_errorbar(aes(ymin=Chla.ng.mm3-se, ymax=Chla.ng.mm3+se), width=.1, position=pd, color="black") + #Error bars
  #  ylim(1,4)+
  xlab("Transplant Site") + ylab(expression("Chlorophyll a " (ng  ~ mm^{-3})))+ #Axis titles
  geom_point(aes(fill=Treatment, shape=Origin), size=14, position=pd, color = "black")+
  scale_shape_manual(values=c(21,24),
                     name = "Reef Zone")+
  scale_fill_manual(values=c("#FFFFFF", "#999999"),
                    name = "Treatment")+ #colour modification
  scale_color_manual(values=c("black", "#999999"),
                     name = "Treatment")+
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/L2018.Chla.vol.pdf", Chla2018Larvae.size, width = 11, height = 11, units = c("in"))


## Statistics

L2018.Chla.meta.size$Origin <- factor(L2018.Chla.meta.size$Origin)
L2018.Chla.meta.size$Treatment <- factor(L2018.Chla.meta.size$Treatment)
L2018.Chla.meta.size$Transplant.Site <- factor(L2018.Chla.meta.size$Transplant.Site)

Chla2018Larvae.2018.anova2 <- lm(Chla.ng.mm3~Origin*Treatment*Transplant.Site, data = L2018.Chla.meta.size)
qqnorm(resid(Chla2018Larvae.2018.anova2))
qqline(resid(Chla2018Larvae.2018.anova2)) 


boxplot(resid(Chla2018Larvae.2018.anova2)~L2018.Chla.meta.size$Origin)
boxplot(resid(Chla2018Larvae.2018.anova2)~L2018.Chla.meta.size$Treatment) 
boxplot(resid(Chla2018Larvae.2018.anova2)~L2018.Chla.meta.size$Transplant.Site)

anova(Chla2018Larvae.2018.anova2)

capture.output(anova(Chla2018Larvae.2018.anova2), file = "output/Statistics/L2018.Chla.vol.csv")



