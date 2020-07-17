#Title: Thermal Transplant 2017-2018 - Total Protein (Adult and Larvae)
#Author: KH Wong
#Date Last Modified: 20200716
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

### Adult 2017 Total Protein ###

# Import Data
A2017.TP.Run1 <- read.csv("data/2017/Protein/Raw_files/TProtein_Adult_20190122_Run1.csv")
A2017.TP.Run2 <- read.csv("data/2017/Protein/Raw_files/TProtein_Adult_20190122_Run2.csv")
A2017.TP.Run3 <- read.csv("data/2017/Protein/Raw_files/TProtein_Adult_20190123_Run3.csv")
A2017.TP.Run4 <- read.csv("data/2017/Protein/Raw_files/TProtein_Adult_20190123_Run4.csv")
A2017.TP.Run5 <- read.csv("data/2017/Protein/Raw_files/TProtein_Adult_20190123_Run5.csv")
A2017.TP.Run6 <- read.csv("data/2017/Protein/Raw_files/TProtein_Adult_20190124_Run6.csv")
A2017.TP.Run7 <- read.csv("data/2017/Protein/Raw_files/TProtein_Adult_20190124_Run7.csv")
A2017.TP.Run8 <- read.csv("data/2017/Protein/Raw_files/TProtein_Adult_20190125_Run8.csv")
A2017.TP.Run9 <- read.csv("data/2017/Protein/Raw_files/TProtein_Adult_20190125_Run9.csv")
A2017.TP.Run10 <- read.csv("data/2017/Protein/Raw_files/TProtein_Adult_20190125_Run10.csv")

A2017.vial.meta <- read.csv("data/2017/Metadata/BIOS2017_Adult_Vial.csv") #Vial metadata
A2017.frag <- read.csv("data/2017/Surface.Area/Adult.Frag.2017.Calculated.csv") #Fragment metadata
A2017.meta <- read.csv("data/2017/Metadata/Colony.Metadata.csv") #Metadata
A2017.TP.well <- read.csv("data/2017/Protein/Adult_Protein_Well_Meta.csv")

# Adding Run numbers into datasets
A2017.TP.Run1$Run <- 1
A2017.TP.Run2$Run <- 2
A2017.TP.Run3$Run <- 3
A2017.TP.Run4$Run <- 4
A2017.TP.Run5$Run <- 5
A2017.TP.Run6$Run <- 6
A2017.TP.Run7$Run <- 7
A2017.TP.Run8$Run <- 8
A2017.TP.Run9$Run <- 9
A2017.TP.Run10$Run <- 10

# Make a unique column for merging
A2017.TP.Run1$run.well <- paste(A2017.TP.Run1$Run, A2017.TP.Run1$Well, sep = "-")
A2017.TP.Run2$run.well <- paste(A2017.TP.Run2$Run, A2017.TP.Run2$Well, sep = "-")
A2017.TP.Run3$run.well <- paste(A2017.TP.Run3$Run, A2017.TP.Run3$Well, sep = "-")
A2017.TP.Run4$run.well <- paste(A2017.TP.Run4$Run, A2017.TP.Run4$Well, sep = "-")
A2017.TP.Run5$run.well <- paste(A2017.TP.Run5$Run, A2017.TP.Run5$Well, sep = "-")
A2017.TP.Run6$run.well <- paste(A2017.TP.Run6$Run, A2017.TP.Run6$Well, sep = "-")
A2017.TP.Run7$run.well <- paste(A2017.TP.Run7$Run, A2017.TP.Run7$Well, sep = "-")
A2017.TP.Run8$run.well <- paste(A2017.TP.Run8$Run, A2017.TP.Run8$Well, sep = "-")
A2017.TP.Run9$run.well <- paste(A2017.TP.Run9$Run, A2017.TP.Run9$Well, sep = "-")
A2017.TP.Run10$run.well <- paste(A2017.TP.Run10$Run, A2017.TP.Run10$Well, sep = "-")

A2017.TP.well$run.well <- paste(A2017.TP.well$Run, A2017.TP.well$Well, sep = "-")

# Merge with metadata
A2017.TP.Run1 <- merge(A2017.TP.Run1, A2017.TP.well, by = "run.well")
A2017.TP.Run2 <- merge(A2017.TP.Run2, A2017.TP.well, by = "run.well")
A2017.TP.Run3 <- merge(A2017.TP.Run3, A2017.TP.well, by = "run.well")
A2017.TP.Run4 <- merge(A2017.TP.Run4, A2017.TP.well, by = "run.well")
A2017.TP.Run5 <- merge(A2017.TP.Run5, A2017.TP.well, by = "run.well")
A2017.TP.Run6 <- merge(A2017.TP.Run6, A2017.TP.well, by = "run.well")
A2017.TP.Run7 <- merge(A2017.TP.Run7, A2017.TP.well, by = "run.well")
A2017.TP.Run8 <- merge(A2017.TP.Run8, A2017.TP.well, by = "run.well")
A2017.TP.Run9 <- merge(A2017.TP.Run9, A2017.TP.well, by = "run.well")
A2017.TP.Run10 <- merge(A2017.TP.Run10, A2017.TP.well, by = "run.well")

# Subtract blanks means for each run (Runs 1, 3, 6, 8 have standards and blanks)
A2017.TP.standardblank1 <- A2017.TP.Run1 %>% 
  filter(Sample.Type == "Blank") %>%
  summarise(blk.avg = mean(X562))

A2017.TP.standardblank3 <- A2017.TP.Run3 %>% 
  filter(Sample.Type == "Blank") %>%
  summarise(blk.avg = mean(X562))

A2017.TP.standardblank6 <- A2017.TP.Run6 %>% 
  filter(Sample.Type == "Blank") %>%
  summarise(blk.avg = mean(X562))

A2017.TP.standardblank8 <- A2017.TP.Run8 %>% 
  filter(Sample.Type == "Blank") %>%
  summarise(blk.avg = mean(X562))

A2017.TP.Run1$abs.corr <- A2017.TP.Run1$X562 - A2017.TP.standardblank1$blk.avg
A2017.TP.Run2$abs.corr <- A2017.TP.Run2$X562 - A2017.TP.standardblank1$blk.avg

A2017.TP.Run3$abs.corr <- A2017.TP.Run3$X562 - A2017.TP.standardblank3$blk.avg
A2017.TP.Run4$abs.corr <- A2017.TP.Run4$X562 - A2017.TP.standardblank3$blk.avg
A2017.TP.Run5$abs.corr <- A2017.TP.Run5$X562 - A2017.TP.standardblank3$blk.avg

A2017.TP.Run6$abs.corr <- A2017.TP.Run6$X562 - A2017.TP.standardblank6$blk.avg
A2017.TP.Run7$abs.corr <- A2017.TP.Run7$X562 - A2017.TP.standardblank6$blk.avg

A2017.TP.Run8$abs.corr <- A2017.TP.Run8$X562 - A2017.TP.standardblank8$blk.avg
A2017.TP.Run9$abs.corr <- A2017.TP.Run9$X562 - A2017.TP.standardblank8$blk.avg
A2017.TP.Run10$abs.corr <- A2017.TP.Run10$X562 - A2017.TP.standardblank8$blk.avg

# Run standards
A2017.TP.standard1 <- A2017.TP.Run1 %>% 
  filter(Sample.Type == "Standard") 

A2017.TP.plot.S1<- ggplot(data = A2017.TP.standard1, aes(x=Concentration, y=abs.corr))+
  ylab("Absorbance (nm)")+ xlab("Concentration") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

A2017.TP.lmstandard1 <- lm (Concentration ~ abs.corr, data = A2017.TP.standard1)
A2017.TP.lmsummary1 <- summary(A2017.TP.lmstandard1)


A2017.TP.standard3 <- A2017.TP.Run3 %>% 
  filter(Sample.Type == "Standard") 

A2017.TP.plot.S3<- ggplot(data = A2017.TP.standard3, aes(x=Concentration, y=abs.corr))+
  ylab("Absorbance (nm)")+ xlab("Concentration") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

A2017.TP.lmstandard3 <- lm (Concentration ~ abs.corr, data = A2017.TP.standard3)
A2017.TP.lmsummary3 <- summary(A2017.TP.lmstandard3)


A2017.TP.standard6 <- A2017.TP.Run6 %>% 
  filter(Sample.Type == "Standard") 

A2017.TP.plot.S6<- ggplot(data = A2017.TP.standard6, aes(x=Concentration, y=abs.corr))+
  ylab("Absorbance (nm)")+ xlab("Concentration") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

A2017.TP.lmstandard6 <- lm (Concentration ~ abs.corr, data = A2017.TP.standard6)
A2017.TP.lmsummary6 <- summary(A2017.TP.lmstandard6)


A2017.TP.standard8 <- A2017.TP.Run8 %>% 
  filter(Sample.Type == "Standard") 

A2017.TP.plot.S8<- ggplot(data = A2017.TP.standard8, aes(x=Concentration, y=abs.corr))+
  ylab("Absorbance (nm)")+ xlab("Concentration") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

A2017.TP.lmstandard8 <- lm (Concentration ~ abs.corr, data = A2017.TP.standard8)
A2017.TP.lmsummary8 <- summary(A2017.TP.lmstandard8)

# Extrapolate concentration values

A2017.TP.Sample1 <- A2017.TP.Run1 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
A2017.TP.Sample1$ConcentrationS <- predict(A2017.TP.lmstandard1, newdata = A2017.TP.Sample1) #using model to get concentration

A2017.TP.Sample2 <- A2017.TP.Run2 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
A2017.TP.Sample2$ConcentrationS <- predict(A2017.TP.lmstandard1, newdata = A2017.TP.Sample2) #using model to get concentration

A2017.TP.Sample3 <- A2017.TP.Run3 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
A2017.TP.Sample3$ConcentrationS <- predict(A2017.TP.lmstandard3, newdata = A2017.TP.Sample3) #using model to get concentration

A2017.TP.Sample4 <- A2017.TP.Run4 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
A2017.TP.Sample4$ConcentrationS <- predict(A2017.TP.lmstandard3, newdata = A2017.TP.Sample4) #using model to get concentration

A2017.TP.Sample5 <- A2017.TP.Run5 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
A2017.TP.Sample5$ConcentrationS <- predict(A2017.TP.lmstandard3, newdata = A2017.TP.Sample5) #using model to get concentration

A2017.TP.Sample6 <- A2017.TP.Run6 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
A2017.TP.Sample6$ConcentrationS <- predict(A2017.TP.lmstandard6, newdata = A2017.TP.Sample6) #using model to get concentration

A2017.TP.Sample7 <- A2017.TP.Run7 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
A2017.TP.Sample7$ConcentrationS <- predict(A2017.TP.lmstandard6, newdata = A2017.TP.Sample7) #using model to get concentration

A2017.TP.Sample8 <- A2017.TP.Run8 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
A2017.TP.Sample8$ConcentrationS <- predict(A2017.TP.lmstandard8, newdata = A2017.TP.Sample8) #using model to get concentration

A2017.TP.Sample9 <- A2017.TP.Run9 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
A2017.TP.Sample9$ConcentrationS <- predict(A2017.TP.lmstandard8, newdata = A2017.TP.Sample9) #using model to get concentration

A2017.TP.Sample10 <- A2017.TP.Run10 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
A2017.TP.Sample10$ConcentrationS <- predict(A2017.TP.lmstandard8, newdata = A2017.TP.Sample10) #using model to get concentration

# Combining datasets

A2017.TP.all <- rbind(A2017.TP.Sample1, A2017.TP.Sample2, A2017.TP.Sample3, A2017.TP.Sample4, A2017.TP.Sample5, A2017.TP.Sample6, A2017.TP.Sample7, A2017.TP.Sample8, A2017.TP.Sample9, A2017.TP.Sample10)
A2017.TP.all <- A2017.TP.all %>%
  select("Vial", "ConcentrationS")

# Adding metadata

A2017.TP.vial <- merge(A2017.TP.all, A2017.vial.meta, by = "Vial")
A2017.TP.vial <- A2017.TP.vial %>%
  select(Vial, ConcentrationS, coral.id, timepoint)

A2017.TP.vial$coral.tp <- paste(A2017.TP.vial$coral.id, A2017.TP.vial$timepoint)

A2017.frag$coral.tp <- paste(A2017.frag$coral.id, A2017.frag$Timepoint)
A2017.frag <- A2017.frag %>%
  select(coral.tp, Surface.Area, homo.vol)

A2017.TP.meta <- merge(A2017.TP.vial, A2017.frag, by = "coral.tp")

# Standardization 

A2017.TP.meta$Total.Protein.ugcm2 <- (A2017.TP.meta$ConcentrationS * A2017.TP.meta$homo.vol * 1.58)/A2017.TP.meta$Surface.Area #Calculating concentration. 1.58 = Dilution factor of acid+base in sample
A2017.TP.meta$Total.Protein.mgcm2 <- A2017.TP.meta$Total.Protein.ugcm2/1000 

# Merging reefzone data
A2017.TP.final <- merge(A2017.TP.meta, A2017.meta, by ="coral.id")

# Removing duplicate run
A2017.TP.final <- A2017.TP.final %>% #removing vial 700, accidentially ran this extra sample
  filter(Vial != "700")

# Summarizing 
mean.A2017.TP <- summarySE(A2017.TP.final, measurevar="Total.Protein.mgcm2", groupvars=c("treatment", "reef.zone", "timepoint"))
mean.A2017.TP$reef.treatment <- paste(mean.A2017.TP$reef.zone, mean.A2017.TP$treatment)

#Reordering factors
mean.A2017.TP$timepoint <- factor(mean.A2017.TP$timepoint, levels = c("Pre-Treatment","Post-Treatment"))

# T1 base plot 
mean.A2017.col.TP <- summarySE(A2017.TP.final, measurevar="Total.Protein.mgcm2", groupvars=c("coral.id", "treatment", "reef.zone", "timepoint"))
mean.A2017.TP.T1 <- summarySE(mean.A2017.col.TP, measurevar="Total.Protein.mgcm2", groupvars=c("timepoint", "reef.zone"))
mean.A2017.TP.T1$timepoint <- factor(mean.A2017.TP.T1$timepoint, levels = c("Pre-Treatment","Post-Treatment"))

# Plotting
pd <- position_dodge(0.3) # moves object .05 to the left and right
#legend.title <- "Treatment"
TP2017Adult.T1 <- ggplot(mean.A2017.TP.T1, aes(x=timepoint, y=Total.Protein.mgcm2, group=reef.zone)) + 
  #  geom_line(position=pd, color="black", aes(linetype=reef.zone), size = 2)+
  geom_errorbar(aes(ymin=Total.Protein.mgcm2-se, ymax=Total.Protein.mgcm2+se), width=.1, size = 1, position=pd, color="black") + #Error bars
  #  ylim(0,5.0)+
  xlab("Timepoint") + ylab(expression("Total Protein " (mg~cm^{-2}))) + #Axis titles
  geom_point(aes(shape=reef.zone), fill = "black", size=14, position=pd, color = "black")+
  scale_shape_manual(values=c(21,24),
                     name = "Reef Zone")+
  #  scale_fill_manual(values=c("#FFFFFF", "#999999"),
  #                    name = "Treatment")+ #colour modification
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  annotate("rect", xmin=0, xmax=1.5, ymin = 1, ymax=3, alpha = 0.2) + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/A2017.TP.T1.pdf", TP2017Adult.T1, width = 11, height = 11, units = c("in"))

# T1 Statistics

A2017.TP.T1 <- mean.A2017.col.TP %>%
  filter(timepoint == "Pre-Treatment")

A2017.TP.T1.anova <- lm(Total.Protein.mgcm2~reef.zone, data = A2017.TP.T1)
qqnorm(resid(A2017.TP.T1.anova))
qqline(resid(A2017.TP.T1.anova))

boxplot(resid(A2017.TP.T1.anova)~A2017.TP.T1$reef.zone)

t.test(Total.Protein.mgcm2~reef.zone, data = A2017.TP.T1)

capture.output(t.test(Total.Protein.mgcm2~reef.zone, data = A2017.TP.T1), file = "output/Statistics/A2017.TP.T1.csv")


# Plotting T2 
A2017.TP.T2 <- mean.A2017.col.TP %>%
  filter(timepoint == "Post-Treatment")

mean.A2017.TP.Col.T2 <- summarySE(A2017.TP.T2, measurevar="Total.Protein.mgcm2", groupvars=c("coral.id","treatment", "reef.zone"))

# Gross photosynthesis 
mean.A2017.TP.T2 <- summarySE(A2017.TP.T2, measurevar="Total.Protein.mgcm2", groupvars=c("treatment", "reef.zone"))

# Plotting
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Treatment"
TP2017Adult.T2 <- ggplot(mean.A2017.TP.T2, aes(x=treatment, y=Total.Protein.mgcm2, group=reef.zone)) + 
  geom_line(position=pd, color="black", aes(linetype=reef.zone), size = 2)+
  geom_errorbar(aes(ymin=Total.Protein.mgcm2-se, ymax=Total.Protein.mgcm2+se), width=.1, size = 1, position=pd, color="black") + #Error bars
#  ylim(0,1)+
  xlab("Timepoint") + ylab(expression("Total Protein " (mg~cm^{-2})))+ #Axis titles
  geom_point(aes(fill=treatment, shape=reef.zone), size=14, position=pd, color = "black")+
  scale_shape_manual(values=c(21,24),
                     name = "Reef Zone")+
  scale_fill_manual(values=c("#FFFFFF", "#999999"),
                    name = "Treatment")+ #colour modification
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  annotate("rect", xmin=0, xmax=0.75, ymin = 1, ymax=3, alpha = 0.2) + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/A2017.TP.T2.pdf", TP2017Adult.T2, width = 11, height = 11, units = c("in"))


##Statistics

A2017.TP.anova.T2 <- lm(Total.Protein.mgcm2~reef.zone*treatment, data = mean.A2017.TP.Col.T2)
qqnorm(resid(A2017.TP.anova.T2))
qqline(resid(A2017.TP.anova.T2))

boxplot(resid(A2017.TP.anova.T2)~mean.A2017.TP.Col.T2$reef.zone)
boxplot(resid(A2017.TP.anova.T2)~mean.A2017.TP.Col.T2$treatment)

anova(A2017.TP.anova.T2)

capture.output(anova(A2017.TP.anova.T2), file = "output/Statistics/A2017.TP.T2.csv")


### 2017 Larval Total Protein ###

L2017.TP <- read.csv("data/2017/Protein/Raw_files/Total_Protein_Larvae.csv")
L2017.Vial <- read.csv("data/2017/Metadata/Vial.metadata.csv")
A2017.meta <- read.csv("data/2017/Metadata/Colony.Metadata.csv") #Metadata
## Making a new data frame for each run because each run has different standards 
# Subsetting Runs
L2017.TP.Run1 <- L2017.TP %>% 
  filter(Run == "1")
L2017.TP.Run2 <- L2017.TP %>% 
  filter(Run == "2")
L2017.TP.Run3 <- L2017.TP %>%
  filter(Run == "3")
L2017.TP.Run4 <- L2017.TP %>%
  filter(Run == "4")
L2017.TP.Run5 <- L2017.TP %>% 
  filter(Run == "5")

#### Subset standards and making a linear model ####
L2017.TP.standard1 <- L2017.TP.Run1 %>% 
  filter(Sample.Type == "Standard") 

L2017.TP.plot.S1<- ggplot(data = L2017.TP.standard1, aes(x=Concentration, y=abs))+
  ylab("Absorbance (nm)")+ xlab("Concentration") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

L2017.TP.lmstandard1 <- lm (Concentration ~ abs, data = L2017.TP.standard1)
L2017.TP.lmsummary1 <- summary(L2017.TP.lmstandard1)
L2017.TP.lmsummary1 

L2017.TP.standard2 <- L2017.TP.Run2 %>% 
  filter(Sample.Type == "Standard") 

L2017.TP.plot.S2<- ggplot(data = L2017.TP.standard2, aes(x=Concentration, y=abs))+
  ylab("Absorbance (nm)")+ xlab("Concentration") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

L2017.TP.lmstandard2 <- lm (Concentration ~ abs, data = L2017.TP.standard2)
L2017.TP.lmsummary2 <- summary(L2017.TP.lmstandard2)
L2017.TP.lmsummary2

#Messed up standard #3, use standard #2 for #3 run samples

L2017.TP.standard4 <- L2017.TP.Run4 %>% 
  filter(Sample.Type == "Standard") 

L2017.TP.plot.S4<- ggplot(data = L2017.TP.standard4, aes(x=Concentration, y=abs))+
  ylab("Absorbance (nm)")+ xlab("Concentration") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

L2017.TP.lmstandard4 <- lm (Concentration ~ abs, data = L2017.TP.standard4)
L2017.TP.lmsummary4 <- summary(L2017.TP.lmstandard4)
L2017.TP.lmsummary4

L2017.TP.standard5 <- L2017.TP.Run5 %>% 
  filter(Sample.Type == "Standard") 

L2017.TP.plot.S5<- ggplot(data = L2017.TP.standard5, aes(x=Concentration, y=abs))+
  ylab("Absorbance (nm)")+ xlab("Concentration") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

L2017.TP.lmstandard5 <- lm (Concentration ~ abs, data = L2017.TP.standard5)
L2017.TP.lmsummary5 <- summary(L2017.TP.lmstandard4)
L2017.TP.lmsummary5

# Calculating protein concentration for samples
L2017.TP.Sample1 <- L2017.TP.Run1 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
L2017.TP.Sample1$ConcentrationS <- predict(L2017.TP.lmstandard1, newdata = L2017.TP.Sample1) #using model to get concentration
L2017.TP.Sample1$Conc.calcS <- (L2017.TP.Sample1$ConcentrationS * L2017.TP.Sample1$DL1 * (1/20)) #Calculating concentration

L2017.TP.Sample2 <- L2017.TP.Run2 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
L2017.TP.Sample2$ConcentrationS <- predict(L2017.TP.lmstandard2, newdata = L2017.TP.Sample2) #using model to get concentration
L2017.TP.Sample2$Conc.calcS <- (L2017.TP.Sample2$ConcentrationS * L2017.TP.Sample2$DL1 * (1/20)) #Calculating concentration

L2017.TP.Sample3 <- L2017.TP.Run3 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
L2017.TP.Sample3$ConcentrationS <- predict(L2017.TP.lmstandard2, newdata = L2017.TP.Sample3) #using model to get concentration
L2017.TP.Sample3$Conc.calcS <- (L2017.TP.Sample3$ConcentrationS * L2017.TP.Sample3$DL1 * (1/20)) #Calculating concentration

L2017.TP.Sample4 <- L2017.TP.Run4 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
L2017.TP.Sample4$ConcentrationS <- predict(L2017.TP.lmstandard4, newdata = L2017.TP.Sample4) #using model to get concentration
L2017.TP.Sample4$Conc.calcS <- (L2017.TP.Sample4$ConcentrationS * L2017.TP.Sample4$DL1 * (1/20)) #Calculating concentration

L2017.TP.Sample5 <- L2017.TP.Run5 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
L2017.TP.Sample5$ConcentrationS <- predict(L2017.TP.lmstandard5, newdata = L2017.TP.Sample5) #using model to get concentration
L2017.TP.Sample5$Conc.calcS <- (L2017.TP.Sample5$ConcentrationS * L2017.TP.Sample5$DL1 * (1/20)) #Calculating concentration

# Combining all datasets together 
L2017.TP.all <- rbind(L2017.TP.Sample1, L2017.TP.Sample2, L2017.TP.Sample3, L2017.TP.Sample4, L2017.TP.Sample5)

L2017.TP.all.clean <- L2017.TP.all %>% select(c(Vial, Conc.calcS))

L2017.TP.vial <- summarySE(L2017.TP.all.clean, measurevar="Conc.calcS", groupvars="Vial")
L2017.TP.meta <- merge(L2017.TP.vial, L2017.Vial, by= "Vial")
L2017.TP.meta2 <- merge(L2017.TP.meta, A2017.meta, by= "coral.id")

L2017.TP.Patch <- subset(L2017.TP.meta2, reef.zone=="Patch")

#### ug protein per mm3 ####

L2017.TP.Patch$Date.coral.ID <- paste(L2017.TP.Patch$Sample.date, L2017.TP.Patch$coral.id, sep = "-")

L.Size.2017.Patch <- read.csv("data/2017/Larval.Size/Patch_mean_Larvalsize.csv")

L2017.TP.Patch.meta <- merge(L.Size.2017.Patch, L2017.TP.Patch, by = "Date.coral.ID")

L2017.TP.Patch.meta$Conc.calcS.size <- L2017.TP.Patch.meta$Conc.calcS/L2017.TP.Patch.meta$Volume_1
L2017.TP.Patch.meta$Conc.calcS.mg.mm3 <- L2017.TP.Patch.meta$Conc.calcS.size/100

L2017.TP.Patch.meta.1 <- L2017.TP.Patch.meta[-c(43), ] #outlier removal
L2017.TP.Final <- summarySE(L2017.TP.Patch.meta.1, measurevar="Conc.calcS.mg.mm3", groupvars=c("treatment", "reef.zone", "Date.coral.ID"))

#Boxplot with jitter
Box.L2017.TP <- ggplot(L2017.TP.Final, aes(x=treatment, y=Conc.calcS.mg.mm3, fill = treatment)) +
  geom_boxplot(width=.3, outlier.colour=NA, position = position_dodge(width = 0.9)) +
  geom_jitter(position = position_jitter(width = 0.1), size = 4) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  #  stat_summary(fun.y=median, geom="line", position = position_dodge(width = 0.9), aes(group=Parental.Treatment))  + 
  #  stat_summary(fun.y=median, geom="point", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("#FFFFFF", "#999999")) +
  xlab("Parental Treatment") + ylab(expression("Total Protein " (mg ~ mm^{-3}))) + #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black")) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/L2017.TP.Patch.size.box.pdf", Box.L2017.TP, width = 11, height = 11, units = c("in"))

# Statistics
L2017.TP.Patch.anova <- lm(Conc.calcS.mg.mm3~treatment, data = L2017.TP.Final)
qqnorm(resid(L2017.TP.Patch.anova))
qqline(resid(L2017.TP.Patch.anova))

boxplot(resid(L2017.TP.Patch.anova)~L2017.TP.Final$treatment)

t.test(Conc.calcS.mg.mm3~treatment, data = L2017.TP.Final)

capture.output(t.test(Conc.calcS.mg.mm3~treatment, data = L2017.TP.Final), file = "output/Statistics/L2017.TP.Patch.csv")


### Total protein Adult 2018 ###

# Import Data
A2018.TP.Run1 <- read.csv("data/2018/Protein/Raw_files/TProtein_Adult2018_20190426_Run1.csv")
A2018.TP.Run2 <- read.csv("data/2018/Protein/Raw_files/TProtein_Adult2018_20190426_Run2.csv")
A2018.TP.Run3 <- read.csv("data/2018/Protein/Raw_files/TProtein_Adult2018_20190514_Run3.csv")
A2018.TP.Run4 <- read.csv("data/2018/Protein/Raw_files/TProtein_Adult2018_20190514_Run4.csv")
A2018.TP.Run5 <- read.csv("data/2018/Protein/Raw_files/TProtein_Adult2018_20190514_Run5.csv")

A2018.vial.meta <- read.csv("data/2018/Metadata/BIOS2018_Adult_Vial.csv") #Vial metadata
A2018.frag.vol <- read.csv("data/2018/Metadata/BIOS2018_Frag_Vol.csv") #Fragment metadata
A2018.frag.sa <- read.csv("data/2018/Surface.Area/Adult.Frag.2018.Calculated.csv") 
A2018.meta <- read.csv("data/2018/Metadata/Sample_Info_Transp.csv") #Metadata
A2018.TP.well <- read.csv("data/2018/Protein/Adult_Protein_2018_Well.csv")

# Adding Run numbers into datasets
A2018.TP.Run1$Run <- 1
A2018.TP.Run2$Run <- 2
A2018.TP.Run3$Run <- 3
A2018.TP.Run4$Run <- 4
A2018.TP.Run5$Run <- 5

# Make a unique column for merging
A2018.TP.Run1$run.well <- paste(A2018.TP.Run1$Run, A2018.TP.Run1$Well, sep = "-")
A2018.TP.Run2$run.well <- paste(A2018.TP.Run2$Run, A2018.TP.Run2$Well, sep = "-")
A2018.TP.Run3$run.well <- paste(A2018.TP.Run3$Run, A2018.TP.Run3$Well, sep = "-")
A2018.TP.Run4$run.well <- paste(A2018.TP.Run4$Run, A2018.TP.Run4$Well, sep = "-")
A2018.TP.Run5$run.well <- paste(A2018.TP.Run5$Run, A2018.TP.Run5$Well, sep = "-")

A2018.TP.well$run.well <- paste(A2018.TP.well$Run, A2018.TP.well$Well, sep = "-")

# Merge with metadata
A2018.TP.Run1 <- merge(A2018.TP.Run1, A2018.TP.well, by = "run.well")
A2018.TP.Run2 <- merge(A2018.TP.Run2, A2018.TP.well, by = "run.well")
A2018.TP.Run3 <- merge(A2018.TP.Run3, A2018.TP.well, by = "run.well")
A2018.TP.Run4 <- merge(A2018.TP.Run4, A2018.TP.well, by = "run.well")
A2018.TP.Run5 <- merge(A2018.TP.Run5, A2018.TP.well, by = "run.well")

# Subtract blanks means for each run (Runs 1, 3have standards and blanks)
A2018.TP.standardblank1 <- A2018.TP.Run1 %>% 
  filter(Sample.Type == "Blank") %>%
  summarise(blk.avg = mean(X562))

A2018.TP.standardblank3 <- A2018.TP.Run3 %>% 
  filter(Sample.Type == "Blank") %>%
  summarise(blk.avg = mean(X562))

A2018.TP.Run1$abs.corr <- A2018.TP.Run1$X562 - A2018.TP.standardblank1$blk.avg
A2018.TP.Run2$abs.corr <- A2018.TP.Run2$X562 - A2018.TP.standardblank1$blk.avg

A2018.TP.Run3$abs.corr <- A2018.TP.Run3$X562 - A2018.TP.standardblank3$blk.avg
A2018.TP.Run4$abs.corr <- A2018.TP.Run4$X562 - A2018.TP.standardblank3$blk.avg
A2018.TP.Run5$abs.corr <- A2018.TP.Run5$X562 - A2018.TP.standardblank3$blk.avg

# Run standards
A2018.TP.standard1 <- A2018.TP.Run1 %>% 
  filter(Sample.Type == "Standard") 
A2018.TP.plot.S1<- ggplot(data = A2018.TP.standard1, aes(x=Concentration, y=abs.corr))+
  ylab("Absorbance (nm)")+ xlab("Concentration") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
A2018.TP.lmstandard1 <- lm (Concentration ~ abs.corr, data = A2018.TP.standard1)
A2018.TP.lmsummary1 <- summary(A2018.TP.lmstandard1)


A2018.TP.standard3 <- A2018.TP.Run3 %>% 
  filter(Sample.Type == "Standard") 
A2018.TP.plot.S3<- ggplot(data = A2018.TP.standard3, aes(x=Concentration, y=abs.corr))+
  ylab("Absorbance (nm)")+ xlab("Concentration") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
A2018.TP.lmstandard3 <- lm (Concentration ~ abs.corr, data = A2018.TP.standard3)
A2018.TP.lmsummary3 <- summary(A2018.TP.lmstandard3)

# Extrapolate concentration values

A2018.TP.Sample1 <- A2018.TP.Run1 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
A2018.TP.Sample1$ConcentrationS <- predict(A2018.TP.lmstandard1, newdata = A2018.TP.Sample1) #using model to get concentration

A2018.TP.Sample2 <- A2018.TP.Run2 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
A2018.TP.Sample2$ConcentrationS <- predict(A2018.TP.lmstandard1, newdata = A2018.TP.Sample2) #using model to get concentration

A2018.TP.Sample3 <- A2018.TP.Run3 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
A2018.TP.Sample3$ConcentrationS <- predict(A2018.TP.lmstandard3, newdata = A2018.TP.Sample3) #using model to get concentration

A2018.TP.Sample4 <- A2018.TP.Run4 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
A2018.TP.Sample4$ConcentrationS <- predict(A2018.TP.lmstandard3, newdata = A2018.TP.Sample4) #using model to get concentration

A2018.TP.Sample5 <- A2018.TP.Run5 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
A2018.TP.Sample5$ConcentrationS <- predict(A2018.TP.lmstandard3, newdata = A2018.TP.Sample5) #using model to get concentration

# Combining datasets

A2018.TP.all <- rbind(A2018.TP.Sample1, 
                      A2018.TP.Sample2, 
                      A2018.TP.Sample3, 
                      A2018.TP.Sample4, 
                      A2018.TP.Sample5)

A2018.TP.all <- A2018.TP.all %>%
  select("Vial", "ConcentrationS")

names(A2018.frag.sa)[2] <- "coral.id"

A2018.TP.coralid <- merge(A2018.TP.all, A2018.vial.meta, by = "Vial")
A2018.TP.vol <- merge(A2018.TP.coralid, A2018.frag.vol, by = "coral.id")
A2018.TP.sa <- merge(A2018.TP.vol, A2018.frag.sa, by = "coral.id")
A2018.TP.meta <- merge(A2018.TP.sa, A2018.meta, by = "coral.id")
A2018.TP.meta <- A2018.TP.meta %>% 
  select(coral.id, ConcentrationS, Homo.Vol.mL, Surface.Area, Origin, Treatment, Transplant.Site)

# Standardization 

A2018.TP.meta$Total.Protein.ugcm2 <- (A2018.TP.meta$ConcentrationS * A2018.TP.meta$Homo.Vol.mL * 1.58)/A2018.TP.meta$Surface.Area #Calculating concentration. 1.58 = Dilution factor of acid+base in sample
A2018.TP.meta$Total.Protein.mgcm2 <- A2018.TP.meta$Total.Protein.ugcm2/1000 

# Summarizing 
mean.TP.Col.A2018 <- summarySE(A2018.TP.meta, measurevar="Total.Protein.mgcm2", groupvars=c("coral.id", "Origin", "Treatment", "Transplant.Site"))
mean.TP.Col.A2018$coral.id <- factor(mean.TP.Col.A2018$coral.id)
mean.TP.Col.A2018$Origin <- factor(mean.TP.Col.A2018$Origin)
mean.TP.Col.A2018$Treatment <- factor(mean.TP.Col.A2018$Treatment)
mean.TP.Col.A2018$Transplant.Site <- factor(mean.TP.Col.A2018$Transplant.Site)

write.csv(mean.TP.Col.A2018, file ="data/2018/Protein/A.Protein.2018.csv")

mean.TP.A2018 <- summarySE(mean.TP.Col.A2018, measurevar="Total.Protein.mgcm2", groupvars=c("Origin", "Treatment", "Transplant.Site"))
mean.TP.A2018$reef.treatment <- paste(mean.TP.A2018$Origin, mean.TP.A2018$Treatment)

#Plotting 
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Treatment"
TP2018Adult <- ggplot(mean.TP.A2018, aes(x=Transplant.Site, y=Total.Protein.mgcm2, group=reef.treatment)) + 
  geom_line(position=pd, aes(linetype=Origin, color = Treatment), size = 2)+
  geom_errorbar(aes(ymin=Total.Protein.mgcm2-se, ymax=Total.Protein.mgcm2+se), width=.1, position=pd, color="black") + #Error bars
  ylim(0.5,2.2)+
  xlab("Transplant Site") + ylab(expression("Total Protein " (mg~cm^{-2}))) + #Axis titles
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

TP2018Adult
ggsave(file = "output/Graphs/A2018.TP.pdf", TP2018Adult,  width = 11, height = 11, units = c("in"))

## Statistics
TP2018Adult.2018.anova <- lm(Total.Protein.mgcm2~Origin*Treatment*Transplant.Site, data = mean.TP.Col.A2018)
qqnorm(resid(TP2018Adult.2018.anova))
qqline(resid(TP2018Adult.2018.anova)) 

boxplot(resid(TP2018Adult.2018.anova)~mean.TP.Col.A2018$Origin)
boxplot(resid(TP2018Adult.2018.anova)~mean.TP.Col.A2018$Treatment) 
boxplot(resid(TP2018Adult.2018.anova)~mean.TP.Col.A2018$Transplant.Site)

anova(TP2018Adult.2018.anova)

capture.output(anova(TP2018Adult.2018.anova), file = "output/Statistics/A2018.TP.csv")



### Larval Total Protein 2018 ### 

# Import Data
L2018.TP.Run1 <- read.csv("data/2018/Protein/Raw_files/TProtein_Larvae2018_20190515_Run1.csv")
L2018.TP.Run2 <- read.csv("data/2018/Protein/Raw_files/TProtein_Larvae2018_20190515_Run2.csv")
L2018.TP.Run3 <- read.csv("data/2018/Protein/Raw_files/TProtein_Larvae2018_20190515_Run3.csv")

L2018.vial.meta <- read.csv("data/2018/Metadata/Larval_2018_vial_metadata.csv") #Vial metadata
L2018.TP.well <- read.csv("data/2018/Protein/Larval_Protein_2018_Well.csv")
A2018.meta <- read.csv("data/2018/Metadata/Sample_Info_Transp.csv")


# Adding Run numbers into datasets
L2018.TP.Run1$Run <- 1
L2018.TP.Run2$Run <- 2
L2018.TP.Run3$Run <- 3

# Make a unique column for merging
L2018.TP.Run1$run.well <- paste(L2018.TP.Run1$Run, L2018.TP.Run1$Well, sep = "-")
L2018.TP.Run2$run.well <- paste(L2018.TP.Run2$Run, L2018.TP.Run2$Well, sep = "-")
L2018.TP.Run3$run.well <- paste(L2018.TP.Run3$Run, L2018.TP.Run3$Well, sep = "-")

L2018.TP.well$run.well <- paste(L2018.TP.well$Run, L2018.TP.well$Well, sep = "-")

# Merge with metadata
L2018.TP.Run1 <- merge(L2018.TP.Run1, L2018.TP.well, by = "run.well")
L2018.TP.Run2 <- merge(L2018.TP.Run2, L2018.TP.well, by = "run.well")
L2018.TP.Run3 <- merge(L2018.TP.Run3, L2018.TP.well, by = "run.well")

# Subtract blanks means for each run (Runs 1, 3have standards and blanks)
L2018.TP.standardblank1 <- L2018.TP.Run1 %>% 
  filter(Sample.Type == "Blank") %>%
  summarise(blk.avg = mean(X562))

L2018.TP.Run1$abs.corr <- L2018.TP.Run1$X562 - L2018.TP.standardblank1$blk.avg
L2018.TP.Run2$abs.corr <- L2018.TP.Run2$X562 - L2018.TP.standardblank1$blk.avg
L2018.TP.Run3$abs.corr <- L2018.TP.Run3$X562 - L2018.TP.standardblank1$blk.avg

# Run standards
L2018.TP.standard1 <- L2018.TP.Run1 %>% 
  filter(Sample.Type == "Standard") 
L2018.TP.plot.S1<- ggplot(data = L2018.TP.standard1, aes(x=Concentration, y=abs.corr))+
  ylab("Absorbance (nm)")+ xlab("Concentration") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
L2018.TP.lmstandard1 <- lm (Concentration ~ abs.corr, data = L2018.TP.standard1)
L2018.TP.lmsummary1 <- summary(L2018.TP.lmstandard1)

# Extrapolate concentration values

L2018.TP.Sample1 <- L2018.TP.Run1 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
L2018.TP.Sample1$ConcentrationS <- predict(L2018.TP.lmstandard1, newdata = L2018.TP.Sample1) #using model to get concentration

L2018.TP.Sample2 <- L2018.TP.Run2 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
L2018.TP.Sample2$ConcentrationS <- predict(L2018.TP.lmstandard1, newdata = L2018.TP.Sample2) #using model to get concentration

L2018.TP.Sample3 <- L2018.TP.Run3 %>% #subsetting Samples
  filter(Sample.Type == "Sample") 
L2018.TP.Sample3$ConcentrationS <- predict(L2018.TP.lmstandard1, newdata = L2018.TP.Sample3) #using model to get concentration

# Combining datasets

L2018.TP.all <- rbind(L2018.TP.Sample1, 
                      L2018.TP.Sample2, 
                      L2018.TP.Sample3)

L2018.TP.all <- L2018.TP.all %>%
  select("Vial", "ConcentrationS")

L2018.TP.coralid <- merge(L2018.TP.all, L2018.vial.meta, by = "Vial")
L2018.TP.meta <- merge(L2018.TP.coralid, A2018.meta, by = "coral.id")
L2018.TP.meta <- L2018.TP.meta %>% 
  select(coral.id, Date.x, ConcentrationS, Sample.Size, Origin, Treatment.y, Transplant.Site)

# Standardization (DL = 4)
L2018.TP.meta$Conc.calcS <- (L2018.TP.meta$ConcentrationS * 4 / L2018.TP.meta$Sample.Size) #Calculating concentration

# Summarizing per larva
mean.TP.col.L2018 <- summarySE(L2018.TP.meta, measurevar="Conc.calcS", groupvars=c("coral.id","Origin", "Treatment.y", "Transplant.Site"))

mean.TP.L2018 <- summarySE(mean.TP.col.L2018, measurevar="Conc.calcS", groupvars=c("Origin", "Treatment.y", "Transplant.Site"))
mean.TP.L2018$reef.treatment <- paste(mean.TP.L2018$Origin, mean.TP.L2018$Treatment)

#### ug protein per mm3 ####

L2018.TP.meta$Date.coral.ID <- paste(L2018.TP.meta$Date.x, L2018.TP.meta$coral.id, sep = "-")

L.Size.2018 <- read.csv("data/2018/Larval.Size/Mean_Larvalsize_ColDay.csv")

L2018.TP.meta.size <- merge(L.Size.2018, L2018.TP.meta, by = "Date.coral.ID")

L2018.TP.meta.size$Conc.calcS.ug.mm3 <- L2018.TP.meta.size$Conc.calcS/L2018.TP.meta.size$Volume
L2018.TP.meta.size$Conc.calcS.mg.mm3 <- L2018.TP.meta.size$Conc.calcS.ug.mm3/100

# Summarizing 
L2018.TP.mean.size.col.day <- summarySE(L2018.TP.meta.size, measurevar="Conc.calcS.mg.mm3", groupvars=c("Date", "coral.id.x", "Origin", "Treatment.y", "Transplant.Site"))

L2018.TP.mean.size <- summarySE(L2018.TP.mean.size.col.day, measurevar="Conc.calcS.mg.mm3", groupvars=c("Origin", "Treatment.y", "Transplant.Site"))

L2018.TP.mean.size$reef.treatment <- paste(L2018.TP.mean.size$Origin, L2018.TP.mean.size$Treatment.y)

#Plotting 
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Treatment"
TP2018Larvae.size <- ggplot(L2018.TP.mean.size, aes(x=Transplant.Site, y=Conc.calcS.mg.mm3, group=reef.treatment)) + 
  geom_line(position=pd, aes(linetype=Origin, color = Treatment.y), size = 2)+
  geom_errorbar(aes(ymin=Conc.calcS.mg.mm3-se, ymax=Conc.calcS.mg.mm3+se), width=.1, position=pd, color="black") + #Error bars
  #  ylim(1,4)+
  xlab("Transplant Site") + ylab(expression(" Total Protein " (mu*g ~ mm^{-3})))+ #Axis titles
  geom_point(aes(fill=Treatment.y, shape=Origin), size=14, position=pd, color = "black")+
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

ggsave(file = "output/Graphs/L2018.TP.vol.pdf", TP2018Larvae.size, width = 11, height = 11, units = c("in"))
