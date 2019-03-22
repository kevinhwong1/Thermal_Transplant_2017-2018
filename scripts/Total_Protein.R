#Title: Thermal Transplant 2017-2018 - Total Protein (Adult and Larvae)
#Author: KH Wong
#Date Last Modified: 20190320
#See Readme file for details 

library(dbplyr)
library(tidyverse)
library(readr)
library(stringr)

### Adult 2017 Total Protein ###

# Import Data
A2017.TP.Run1 <- read.csv("data/2017/Protein/TProtein_Adult_20190122_Run1.csv")
A2017.TP.Run2 <- read.csv("data/2017/Protein/TProtein_Adult_20190122_Run2.csv")
A2017.TP.Run3 <- read.csv("data/2017/Protein/TProtein_Adult_20190123_Run3.csv")
A2017.TP.Run4 <- read.csv("data/2017/Protein/TProtein_Adult_20190123_Run4.csv")
A2017.TP.Run5 <- read.csv("data/2017/Protein/TProtein_Adult_20190123_Run5.csv")
A2017.TP.Run6 <- read.csv("data/2017/Protein/TProtein_Adult_20190124_Run6.csv")
A2017.TP.Run7 <- read.csv("data/2017/Protein/TProtein_Adult_20190124_Run7.csv")
A2017.TP.Run8 <- read.csv("data/2017/Protein/TProtein_Adult_20190125_Run8.csv")
A2017.TP.Run9 <- read.csv("data/2017/Protein/TProtein_Adult_20190125_Run9.csv")
A2017.TP.Run10 <- read.csv("data/2017/Protein/TProtein_Adult_20190125_Run10.csv")

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

# Plotting
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Treatment"
TP2017Adult <- ggplot(mean.A2017.TP, aes(x=timepoint, y=Total.Protein.mgcm2, group=reef.treatment)) + 
  geom_line(position=pd, color="black")+
  geom_errorbar(aes(ymin=Total.Protein.mgcm2-se, ymax=Total.Protein.mgcm2+se), width=.1, position=pd, color="black") + #Error bars
  ylim(1.3,2.6)+
  xlab("Timepoint") + ylab(expression("Total Protein " (mu*g~cm^{-2})))+ #Axis titles
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

ggsave(file = "output/Graphs/A2017.TP.pdf", TP2017Adult)

# Statsitics 
A2017.TP.anova <- lm(Total.Protein.mgcm2~reef.zone*treatment*timepoint, data = A2017.TP.final)
qqnorm(resid(A2017.TP.anova))
qqline(resid(A2017.TP.anova))

boxplot(resid(A2017.TP.anova)~A2017.TP.final$reef.zone)
boxplot(resid(A2017.TP.anova)~A2017.TP.final$treatment)
boxplot(resid(A2017.TP.anova)~A2017.TP.final$timepoint)
anova(A2017.TP.anova)

capture.output(anova(A2017.TP.anova), file = "output/Statistics/A2017.TP.csv")


### 2017 Larval Total Protein ###

L2017.TP <- read.csv("data/2017/Protein/Total_Protein_Larvae.csv")
L2017.Vial <- read.csv("data/2017/Metadata/Vial.metadata.csv")

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

##### ug Protein per larva #####
#Getting means
L2017.TP.Patch.mean <- summarySE(L2017.TP.Patch, measurevar="Conc.calcS", groupvars="treatment")

#Plotting reaction norm
pd <- position_dodge(0.1) # moves object .05 to the left and right
legend.title <- "Treatment"
L2017.TP.Patch.plot<- ggplot(L2017.TP.Patch.mean, aes(x=treatment, y=Conc.calcS)) + 
  geom_line(position=pd, color="black")+
  geom_errorbar(aes(ymin=Conc.calcS-se, ymax=Conc.calcS+se), width=.1, position=pd, color="black") + #Error bars
  ylim(30,45)+
  xlab("Treatment") + ylab(expression(" Total Protein " (mu*g ~ Larva^{-1})))+ #Axis titles
  geom_point(aes(fill=treatment), color ="black", pch=21, size=4)+
  scale_fill_manual(legend.title, values=c("dodgerblue3", "tomato1"))+ #colour modification
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.title.x = element_blank())


L2017.TP.Patch.plot
ggsave(file = "output/Graphs/L2017.TP.Patch.pdf", L2017.TP.Patch.plot)


# Statistics
L2017.TP.Patch.anova <- lm(Conc.calcS~treatment, data = L2017.TP.Patch)
qqnorm(resid(L2017.TP.Patch.anova))
qqline(resid(L2017.TP.Patch.anova))

boxplot(resid(L2017.TP.Patch.anova)~L2017.TP.Patch$treatment)

t.test(Conc.calcS~treatment, data = L2017.TP.Patch)

capture.output(t.test(Conc.calcS~treatment, data = L2017.TP.Patch), file = "output/Statistics/L2017.TP.Patch.csv")
