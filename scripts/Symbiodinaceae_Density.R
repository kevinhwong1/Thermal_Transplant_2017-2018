#Title: Thermal Transplant 2017-2018 - Symbiodinaceae Densitiy (Adult and Larvae)
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

write.csv(A2017.zoox.meta, 'data/2017/Zoox/A2017.zoox.calc.csv')

#Summarizing
A2017.zoox.mean <- summarySE(A2017.zoox.meta, measurevar="Cells.cm2.x6", groupvars=c("treatment", "reef.zone", "timepoint"))
A2017.zoox.mean$reef.treatment <- paste(A2017.zoox.mean$reef.zone, A2017.zoox.mean$treatment)

#Reordering factors
A2017.zoox.mean$timepoint <- factor(A2017.zoox.mean$timepoint, levels = c("Pre-Treatment","Post-Treatment"))

#Plotting T1 (Pre-treatment)
A2017.zoox.mean2 <- summarySE(A2017.zoox.meta, measurevar="Cells.cm2.x6", groupvars=c("reef.zone", "timepoint"))

#Reordering factors
A2017.zoox.mean2$timepoint <- factor(A2017.zoox.mean2$timepoint, levels = c("Pre-Treatment","Post-Treatment"))
pd <- position_dodge(0.3) # moves object .05 to the left and right

Zoox2017Adult.T1 <- ggplot(A2017.zoox.mean2, aes(x=timepoint, y=Cells.cm2.x6, group=reef.zone)) + 
  geom_errorbar(aes(ymin=Cells.cm2.x6-se, ymax=Cells.cm2.x6+se), width=.1, size = 1, position=pd, color="black") + #Error bars
  ylim(0,1.6)+
  xlab("Timepoint") + ylab(expression("Cell Density " (10^{6} ~ cm^{-2})))+ #Axis titles
  geom_point(aes(shape=reef.zone), fill = "black", size=14, position=pd, color = "black")+
  scale_shape_manual(values=c(21,24),
                     name = "Reef Zone")+
  annotate("rect", xmin=0, xmax=1.5, ymin = 0, ymax=1.5, alpha = 0.2) + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/A2017.Zoox.T1.pdf", Zoox2017Adult.T1, width = 11, height = 11, units = c("in"))

# T1 Statistics
A2017.zoox.T1 <- A2017.zoox.meta %>%
  filter(timepoint == "Pre-Treatment")

A2017.zoox.T1.anova <- lm(Cells.cm2.x6~reef.zone, data = A2017.zoox.T1)
qqnorm(resid(A2017.zoox.T1.anova))
qqline(resid(A2017.zoox.T1.anova))

boxplot(resid(A2017.zoox.T1.anova)~A2017.zoox.T1$reef.zone)

t.test(Cells.cm2.x6~reef.zone, data = A2017.zoox.T1)

capture.output(t.test(Cells.cm2.x6~reef.zone, data = A2017.zoox.T1), file = "output/Statistics/A2017.Zoox.T1.csv")


##T2 analysis 

A2017.zoox.T2 <- A2017.zoox.meta %>%
  filter(timepoint == "Post-Treatment")

# Summarizing 
mean.A2017.zoox.T2 <- summarySE(A2017.zoox.T2, measurevar="Cells.cm2.x6", groupvars=c("treatment", "reef.zone"))

# Plotting
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Treatment"
Zoox2017Adult.T2 <- ggplot(mean.A2017.zoox.T2, aes(x=treatment, y=Cells.cm2.x6, group=reef.zone)) + 
  geom_line(position=pd, color="black", aes(linetype=reef.zone), size = 2)+
  geom_errorbar(aes(ymin=Cells.cm2.x6-se, ymax=Cells.cm2.x6+se), width=.1, size = 1, position=pd, color="black") + #Error bars
    ylim(0,1.6)+
  xlab("Timepoint") + ylab(expression("Cell Density " (10^{6} ~ cm^{-2})))+ #Axis titles
  geom_point(aes(fill=treatment, shape=reef.zone), size=14, position=pd, color = "black")+
  scale_shape_manual(values=c(21,24),
                     name = "Reef Zone")+
  scale_fill_manual(values=c("#FFFFFF", "#999999"),
                    name = "Treatment")+ #colour modification
  annotate("rect", xmin=0, xmax=0.75, ymin = 0, ymax=1.5, alpha = 0.2) + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/A2017.Zoox.T2.pdf", Zoox2017Adult.T2, width = 11, height = 11, units = c("in"))

## Statistics
zoox.A2017.anova <- lm(Cells.cm2.x6~reef.zone*treatment, data = A2017.zoox.T2)
qqnorm(resid(zoox.A2017.anova))
qqline(resid(zoox.A2017.anova))

boxplot(resid(zoox.A2017.anova)~A2017.zoox.T2$reef.zone)
boxplot(resid(zoox.A2017.anova)~A2017.zoox.T2$treatment)

anova(zoox.A2017.anova)

capture.output(anova(zoox.A2017.anova), file = "output/Statistics/A2017.Zoox.csv")




### 2017 Larval Symbiodinaceae Densities ###

#Import Data
L2017.zoox <- read.csv("data/2017/Zoox/Larval_CellCounts.csv")
A2017.meta <- read.csv("data/2017/Metadata/Colony.Metadata.csv") #Metadata
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
L2017.zoox.patch.mean <- summarySE(L2017.zoox.patch, measurevar="Cells.larvae.x3", groupvars=c("treatment", "coral.id", "date"))

#Plotting
Box.L2017.Zoox.larva <- ggplot(L2017.zoox.patch.mean, aes(x=treatment, y=Cells.larvae.x3, fill = treatment)) +
  geom_boxplot(width=.3, outlier.colour=NA, position = position_dodge(width = 0.9)) +
  geom_jitter(position = position_jitter(width = 0.1), size =4) +
  #  stat_summary(fun.y=median, geom="line", position = position_dodge(width = 0.9), aes(group=Parental.Treatment))  + 
  #  stat_summary(fun.y=median, geom="point", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("#FFFFFF", "#999999")) +
  xlab("Parental Treatment") + ylab(expression("Cell Density " (10^{3} ~ larva^{-1}))) + #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black")) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/L2017.Zoox.Patch.larva.box.pdf", Box.L2017.Zoox.larva, width = 11, height = 11, units = c("in"))

# Statistics
capture.output(t.test(Cells.larvae.x3~treatment, data = L2017.zoox.patch.mean), file = "output/Statistics/L2017.Zoox.larva.csv")



#### larval zoox per mm3 ####

L2017.zoox.patch$Date.coral.ID <- paste(L2017.zoox.patch$date, L2017.zoox.patch$coral.id, sep = "-")

L.Size.2017.Patch <- read.csv("data/2017/Larval.Size/Patch_mean_Larvalsize.csv")

L2017.zoox.Patch.size <- merge(L.Size.2017.Patch, L2017.zoox.patch, by = "Date.coral.ID")

L2017.zoox.Patch.size$Cells.x3.mm3 <- L2017.zoox.Patch.size$Cells.larvae.x3/L2017.zoox.Patch.size$Volume_1

mean.L2017.zoox.size <- summarySE(L2017.zoox.Patch.size, measurevar="Cells.x3.mm3", groupvars=c("treatment", "reef.zone", "Date.Release", "coral.id.x"))
write.csv(mean.L2017.zoox.size, file = "data/2017/Zoox/L2017_zoox.vol.csv")


#Plotting
Box.L2017.Zoox <- ggplot(mean.L2017.zoox.size, aes(x=treatment, y=Cells.x3.mm3, fill = treatment)) +
  geom_boxplot(width=.3, outlier.colour=NA, position = position_dodge(width = 0.9)) +
  geom_jitter(position = position_jitter(width = 0.1), size =4) +
  #  stat_summary(fun.y=median, geom="line", position = position_dodge(width = 0.9), aes(group=Parental.Treatment))  + 
  #  stat_summary(fun.y=median, geom="point", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("#FFFFFF", "#999999")) +
  xlab("Parental Treatment") + ylab(expression("Cell Density " (10^{3} ~ mm^{-3}))) + #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black")) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/L2017.Zoox.Patch.size.box.pdf", Box.L2017.Zoox, width = 11, height = 11, units = c("in"))

# Statistics
capture.output(t.test(Cells.x3.mm3~treatment, data = mean.L2017.zoox.size), file = "output/Statistics/L2017.Zoox.size.csv")


# Statistics (mixed effect modelling)

### orthogonal random effect design 

#One way model with date and coral ID as random factor
Zoox2017Larvae.lmer <- lmer(Cells.x3.mm3~1+treatment + (1|Date.Release) +  (1|coral.id.x), data = mean.L2017.zoox.size, REML=FALSE)
qqnorm(resid(Zoox2017Larvae.lmer)) # Normality
qqline(resid(Zoox2017Larvae.lmer)) # Normal

boxplot(resid(Zoox2017Larvae.lmer)~mean.L2017.zoox.size$treatment) # Variance 

summary(Zoox2017Larvae.lmer)

#One way model with date as random factor
Zoox2017Larvae.lmer2 <- lmer(Cells.x3.mm3~1+treatment + (1|Date.Release), data = mean.L2017.zoox.size, REML=FALSE)
qqnorm(resid(Zoox2017Larvae.lmer2)) # Normality
qqline(resid(Zoox2017Larvae.lmer2)) # Normal

boxplot(resid(Zoox2017Larvae.lmer2)~mean.L2017.zoox.size$treatment) # Variance 

summary(Zoox2017Larvae.lmer2)

#one way model with coral.id as random factor
Zoox2017Larvae.lmer3 <- lmer(Cells.x3.mm3~1+treatment + (1|coral.id.x), data = mean.L2017.zoox.size, REML=FALSE)
qqnorm(resid(Zoox2017Larvae.lmer3)) # Normality
qqline(resid(Zoox2017Larvae.lmer3)) # Normal

boxplot(resid(Zoox2017Larvae.lmer3)~mean.L2017.zoox.size$treatment) # Variance 

summary(Zoox2017Larvae.lmer3)

## Model Comparisons
lrt(Zoox2017Larvae.lmer2, Zoox2017Larvae.lmer) # Model II is best
lrt(Zoox2017Larvae.lmer3, Zoox2017Larvae.lmer) # Model III is best

#MODEL SELECTION = II
capture.output(Anova(Zoox2017Larvae.lmer2), file = "output/Statistics/L.2017.zoox.lmer.csv")


### 2018 Adult Symbiodinaceae Densities ###

#Import Data
A2018.adult.Zoox <- read.csv("data/2018/Zoox/Adult_Zoox_2018.csv")

A2018.Frag.vol <- read.csv("data/2018/Metadata/BIOS2018_Frag_Vol.csv") #Vial metadata
A2018.Frag.sa <- read.csv("data/2018/Surface.Area/Adult.Frag.2018.Calculated.csv") #Fragment metadata
A2018.meta <- read.csv("data/2018/Metadata/Sample_Info_Transp.csv") #Metadata

#Calculating cells/Larvae from cells/count
A2018.adult.Zoox$Average <- rowMeans(A2018.adult.Zoox[ ,7:12]) #averaging counts
A2018.adult.Zoox$Cells.mL <- (A2018.adult.Zoox$Average/A2018.adult.Zoox$X.squares.count)/0.0001 #number of cells per mL - volume of haemocytometer

#Accounting for counting dilution factor

A2018.adult.Zoox$Cell.ml.DF <- A2018.adult.Zoox$Cells.mL * A2018.adult.Zoox$Dilution.Factor

#Removing metadata columns
A2018.adult.Zoox2 <- A2018.adult.Zoox %>%
  select("coral.id", "Cell.ml.DF")

names(A2018.Frag.sa)[2] <- "coral.id"

#Merging dataframes
A2018.zoox.colony <- merge(A2018.adult.Zoox2, A2018.meta, by = "coral.id")
A2018.zoox.SA <- merge(A2018.zoox.colony, A2018.Frag.sa, by = "coral.id")
A2018.zoox.meta <- merge(A2018.zoox.SA, A2018.Frag.vol, by = "coral.id")

#Removing metadata columns
A2018.zoox.meta <- A2018.zoox.meta %>%
  select("coral.id","Origin", "Treatment", "Transplant.Site", "Cell.ml.DF", "Homo.Vol.mL", "Surface.Area")

#Accounting for homogenate volume 
A2018.zoox.meta$Cells.mL.vol <- A2018.zoox.meta$Cell.ml.DF*A2018.zoox.meta$Homo.Vol.mL 

#Normalizing to fragment surface area
A2018.zoox.meta$Cells.cm2 <- A2018.zoox.meta$Cells.mL.vol/A2018.zoox.meta$Surface.Area 
A2018.zoox.meta$Cells.cm2.x6 <- A2018.zoox.meta$Cells.cm2 / 1000000

write.csv(A2018.zoox.meta, file="data/2018/Zoox/A2018_Zoox_Calc.csv")

#Summarizing
A2018.zoox.meta.mean <- summarySE(A2018.zoox.meta, measurevar="Cells.cm2.x6", groupvars=c("Origin", "Treatment", "Transplant.Site"))


A2018.zoox.meta.mean$reef.treatment <- paste(A2018.zoox.meta.mean$Origin, A2018.zoox.meta.mean$Treatment)

#Plotting reaction norm
pd <- position_dodge(0.1) # moves object .05 to the left and right

zoox.adult.2018<- ggplot(A2018.zoox.meta.mean, aes(x=Transplant.Site, y=Cells.cm2.x6, group=reef.treatment)) + 
  geom_line(position=pd, aes(linetype=Origin, color = Treatment), size = 2)+
  geom_errorbar(aes(ymin=Cells.cm2.x6-se, ymax=Cells.cm2.x6+se), width=.1, position=pd, color="black") + #Error bars
  ylim(0.1,1.6)+
  xlab("Transplant Site") + ylab(expression("Cell Density " (10^{6} ~ cm^{-2})))+ #Axis titles
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

zoox.adult.2018
ggsave(file = "output/Graphs/A2018.Zoox.pdf", zoox.adult.2018, width = 11, height = 11, units = c("in"))


#Statistics
Zoox2018Adult.anova <- aov(Cells.cm2.x6~Origin*Treatment*Transplant.Site, data = A2018.zoox.meta)
qqnorm(resid(Zoox2018Adult.anova))
qqline(resid(Zoox2018Adult.anova)) 

boxplot(resid(Zoox2018Adult.anova)~A2018.zoox.meta$Origin)
boxplot(resid(Zoox2018Adult.anova)~A2018.zoox.meta$Treatment) 
boxplot(resid(Zoox2018Adult.anova)~A2018.zoox.meta$Transplant.Site)

summary(Zoox2018Adult.anova)
anova(Zoox2018Adult.anova)

capture.output(anova(Zoox2018Adult.anova), file = "output/Statistics/A2018.Zoox.csv")

# Post-Hoc 

A.2018.Zoox.PH.OxTrans <- TukeyHSD(Zoox2018Adult.anova, conf.level = 0.95)
capture.output(A.2018.Zoox.PH.OxTrans, file = "output/Statistics/A2018.Zoox.PH.csv")

#PostHoc Tukey adjustment comparison of least-squares means
A.2018.Zoox.OxTrans <- lsmeans(Zoox2018Adult.anova, ~ Origin*Transplant.Site, adjust="tukey") #compute least-squares means for Treatment*Day from ANOVA model
A.2018.Zoox.pairs.OxTrans <- multcomp::cld(A.2018.Zoox.OxTrans, alpha=.05, Letters=letters) #list pairwise tests and letter display
A.2018.Zoox.pairs.OxTrans #view results
capture.output(A.2018.Zoox.pairs.OxTrans, file = "output/Statistics/A.2018.Zoox.pairs.OxTrans.csv")


##### Larval 2018 zoox density #####

#Import Data
L2018.zoox <- read.csv("data/2018/Zoox/Larval_2018_Zoox.csv")

#Cleaning up data frame
L2018.zoox <- L2018.zoox[,colSums(is.na(L2018.zoox))<nrow(L2018.zoox)] #Removing columns with all NAs
L2018.zoox.clean <-L2018.zoox[complete.cases(L2018.zoox), ]

#Calculating cells/Larvae from cells/count
L2018.zoox.clean$Average <- rowMeans(L2018.zoox.clean[ ,11:16]) #averaging counts
L2018.zoox.clean$Cells.mL <- (L2018.zoox.clean$Average/L2018.zoox.clean$X.squares.count)/0.0001 #number of cells per mL - volume of haemocytometer
L2018.zoox.clean$Cells.mL.Diluted <- L2018.zoox.clean$Cells.mL/L2018.zoox.clean$Dilution.Factor #factoring in dilution
L2018.zoox.clean$Cells.larvae <- L2018.zoox.clean$Cells.mL.Diluted/L2018.zoox.clean$Sample.Size #diving by number of larvae per vial

#Remove metadata columns
#L2018.zoox.clean2 <- L2018.zoox.clean %>% select(Date, Colony, Cells.larvae)

#Renaming columns
colnames(L2018.zoox.clean)[colnames(L2018.zoox.clean)=='Colony']<-'coral.id' #changing column name 

A2018.meta <- read.csv("data/2018/Metadata/Sample_Info_Transp.csv")

#Attaching metadata
L2018.zoox.clean3 <- merge(L2018.zoox.clean, A2018.meta, by = "coral.id")

#making x103
L2018.zoox.clean3$Cells.larvae.x3 <-  L2018.zoox.clean3$Cells.larvae / 1000

#Summarizing
L2018.zoox.colday <- summarySE(L2018.zoox.clean3, measurevar="Cells.larvae.x3", groupvars=c("coral.id", "Date.x", "Origin", "Treatment.y", "Transplant.Site"))

L2018.zoox.mean <- summarySE(L2018.zoox.clean3, measurevar="Cells.larvae.x3", groupvars=c("Origin", "Treatment.y", "Transplant.Site"))
L2018.zoox.mean$reef.treatment <- paste(L2018.zoox.mean$Origin, L2018.zoox.mean$Treatment)

#### cellsx3 per mm3 ####

L2018.zoox.clean3$Date.coral.ID <- paste(L2018.zoox.clean3$Date.x, L2018.zoox.clean3$coral.id, sep = "-")

L.Size.2018 <- read.csv("data/2018/Larval.Size/Mean_Larvalsize_ColDay.csv")

L2018.zoox.meta.size <- merge(L.Size.2018, L2018.zoox.clean3, by = "Date.coral.ID")

L2018.zoox.meta.size$Cells.x3.mm3 <- L2018.zoox.meta.size$Cells.larvae.x3/L2018.zoox.meta.size$Volume

# Summarizing 
L2018.zoox.mean.col.size <- summarySE(L2018.zoox.meta.size, measurevar="Cells.x3.mm3", groupvars=c("Date", "coral.id.x", "Origin", "Treatment.y", "Transplant.Site"))
write.csv(L2018.zoox.mean.col.size, file = "data/2018/Zoox/L2018_Zoox_calc.csv")
L2018.zoox.mean.size <- summarySE(L2018.zoox.meta.size, measurevar="Cells.x3.mm3", groupvars=c("Origin", "Treatment.y", "Transplant.Site"))

#mean.TP.L2018 <- summarySE(mean.TP.col.L2018, measurevar="Conc.calcS", groupvars=c("Origin", "Treatment.y", "Transplant.Site"))
L2018.zoox.mean.size$reef.treatment <- paste(L2018.zoox.mean.size$Origin, L2018.zoox.mean.size$Treatment.y)

#Plotting 
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Treatment"
Zoox2018Larvae.size <- ggplot(L2018.zoox.mean.size, aes(x=Transplant.Site, y=Cells.x3.mm3, group=reef.treatment)) + 
  geom_line(position=pd, aes(linetype=Origin, color = Treatment.y), size = 2)+
  geom_errorbar(aes(ymin=Cells.x3.mm3-se, ymax=Cells.x3.mm3+se), width=.1, position=pd, color="black") + #Error bars
  ylim(10,55)+
  xlab("Transplant Site") + ylab(expression("Cell Density " (10^{3} ~ mm^{-3})))+ #Axis titles
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

ggsave(file = "output/Graphs/L2018.Zoox.vol.pdf", Zoox2018Larvae.size, width = 11, height = 11, units = c("in"))


# ## Making residuals
# A.Zoox.2018.TPatch <- A2018.zoox.meta %>%
#   filter(Transplant.Site == "Patch")
# 
# A.Zoox.2018.TRim <- A2018.zoox.meta %>%
#   filter(Transplant.Site == "Rim")
# 
# A.Zoox.2018.TPatch$coral.id <- A.Zoox.2018.TPatch$coral.id %>% str_replace("-A","")
# A.Zoox.2018.TRim$coral.id <- A.Zoox.2018.TRim$coral.id %>% str_replace("-B","")
# 
# A.Zoox.2018.comp <- merge(A.Zoox.2018.TPatch, A.Zoox.2018.TRim, by = "coral.id")
# 
# # residual calculation
# A.Zoox.2018.comp$resid.A.Zoox <- resid(lm(A.Zoox.2018.comp$Cells.cm2.x6.y - A.Zoox.2018.comp$Cells.cm2.x6.x ~0))
# 
# # Summarizing residuals
# A.Zoox.resid.2018.mean <- summarySE(A.Zoox.2018.comp, measurevar="resid.A.Zoox", groupvars=c("Treatment.x", "Origin.x"))
# 
# colnames(A.Zoox.resid.2018.mean) [1] <- "Treatment"
# colnames(A.Zoox.resid.2018.mean) [2] <- "Origin"
# 
# write.csv(A.Zoox.resid.2018.mean, file ="output/Residual_Analysis/Resid.A.Zoox.csv")
# 
# # Statsitics
# A.Zoox.resid.anova <- lm(resid.A.Zoox~Origin.x*Treatment.x, data = A.Zoox.2018.comp)
# qqnorm(resid(A.Zoox.resid.anova))
# qqline(resid(A.Zoox.resid.anova))
# 
# boxplot(resid(A.Zoox.resid.anova)~A.Zoox.2018.comp$Origin.x)
# boxplot(resid(A.Zoox.resid.anova)~A.Zoox.2018.comp$Treatment.x)
# 
# anova(A.Zoox.resid.anova)
# 
# capture.output(anova(A.Zoox.resid.anova), file = "output/Statistics/A2018.Zoox.Resid.csv")
