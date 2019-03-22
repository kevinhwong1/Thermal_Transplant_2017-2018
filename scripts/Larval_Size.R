#Title: Thermal Transplant 2017-2018 - Larval Size
#Author: KH Wong
#Date Last Modified: 20190321
#See Readme file for details 

#Load Libaries
library(dbplyr)
library(tidyverse)
library(readr)

### 2017 Larval Size ###

#load data
L2017.vol.data<-read.csv('data/2017/Larval.Size/TA_Larval_Size.csv', header=T, sep=",") #loading data

colnames(L2017.vol.data)[colnames(L2017.vol.data)=='Colony']<-'coral.id' #changing column name 

L2017.vol.data<- subset(L2017.vol.data, Reef=="Patch")
L2017.vol.T1 <- subset(L2017.vol.data, Timepoint=="T1") #Subsetting T1 from T2 and creating a variable for T1

# Summarizing
L2017.vol.mean <- summarySE(L2017.vol.T1, measurevar="Volume_1", groupvars=c("Treatment.1","Reef")) #summarize vol by treatment and reef

# Plotting
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Parental Treatment"
larval.size.2017 <- ggplot(L2017.vol.mean, aes(x=Treatment.1, y=Volume_1)) + 
  geom_line(position=pd, color="black")+
  geom_errorbar(aes(ymin=Volume_1-se, ymax=Volume_1+se), width=.1, position=pd, color="black") + #Error bars
  ylim(0.2,0.4)+
  xlab("Treatment") + ylab(expression("Larval Volume" ~ (mm^{3})))+ #Axis titles
  geom_point(aes(fill=Treatment.1), color ="black", pch=21, size=4)+
  scale_fill_manual(legend.title, values=c("dodgerblue3", "tomato1"))+ #colour modification
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        axis.title.x = element_blank())+
  theme(legend.position = "none")

ggsave(file = "output/Graphs/L2017.Vol.pdf", larval.size.2017)

# Statistics
L2017.Vol.anova <- lm(Volume_1~Treatment.1, data = L2017.vol.T1)
qqnorm(resid(L2017.Vol.anova))
qqline(resid(L2017.Vol.anova))

boxplot(resid(L2017.Vol.anova)~L2017.vol.T1$Treatment.1)

t.test(Volume_1~Treatment.1, data = L2017.vol.T1)

capture.output(t.test(Volume_1~Treatment.1, data = L2017.vol.T1), file = "output/Statistics/L2017.Vol.csv")

### 2018 Larval Size ###

#load data
L2018.vol.data<-read.csv('data/2018/Larval.Size/Transp_Larval_Size.csv', header=T, sep=",") #loading data
colony.info.2018 <- read.csv('data/2018/Metadata/Sample_Info_Transp.csv')


colnames(L2018.vol.data)[colnames(L2018.vol.data)=='Colony']<-'coral.id' #changing column name 

#Calculating volume
##uequation for volume of an elliptical sphere, V=4/3(pi)ab^2, where a is 1/2 width and b is 1/2 length (Isomura and Nishihira 2001) 

L2018.vol.data$a <- L2018.vol.data$Width/2
L2018.vol.data$b <- L2018.vol.data$Length/2
L2018.vol.data$Volume <- 4/3*pi*L2018.vol.data$a*(L2018.vol.data$b^2) 

#adding sample info
L2018.vol.info <- merge(L2018.vol.data, colony.info.2018, by="coral.id")

#Removing As and Bs
L2018.vol.info$coral.id.orig <- substr(as.character(L2018.vol.info$coral.id), #removing As and Bs in coral.id
                                  start= 1, 
                                  stop= nchar(as.character(L2018.vol.info$coral.id) )-2 ) 

#Unique column
L2018.vol.info$Origin.Treatment <- paste(L2018.vol.info$Origin, L2018.vol.info$Treatment)
L2018.vol.info <- L2018.vol.info %>%
  select(coral.id, Volume, Date.x, Origin, Treatment, Transplant.Site, Group, coral.id.orig, Origin.Treatment)

#Mean Sizes by colony
L2018.vol.mean.col <- summarySE(L2018.vol.info, measurevar="Volume", groupvars=c("coral.id","coral.id.orig","Origin", "Treatment", "Origin.Treatment", "Transplant.Site")) #summarize vol by treatment and reef

#reformat dataframe
L2018.vol.TPatch <- L2018.vol.mean.col %>%
  filter(Transplant.Site == "Patch")

L2018.vol.TRim <- L2018.vol.mean.col %>%
  filter(Transplant.Site == "Rim")

L2018.vol.comp <- merge(L2018.vol.TPatch, L2018.vol.TRim, by = "coral.id.orig")

#Making residuals
L2018.vol.comp$resid.1 <- resid(lm(L2018.vol.comp$Volume.y - L2018.vol.comp$Volume.x ~0))

#Summarizing residuals
L2018.vol.resid.mean <- summarySE(L2018.vol.comp, measurevar="resid.1", groupvars=c("Treatment.x", "Origin.x"))

#Plotting
pd <- position_dodge(0.1) # moves object .05 to the left and right
size.resid.2018<- ggplot(L2018.vol.resid.mean, aes(x=Treatment.x, y=resid.1, group=Origin.x, shape = Origin.x)) + 
  geom_line(position=pd, color="black")+
  geom_errorbar(aes(ymin=resid.1-se, ymax=resid.1+se), width=.1, position=pd, color="black") + #Error bars
  geom_hline(yintercept = 0,linetype="dashed") +
  #  ylim(0,16)+
  xlab("Treatment") + ylab(expression(paste("Larval Volume ", cm^{3}, " ( Rim - Patch)"))) + #Axis titles
  geom_point(aes(shape=Origin.x), position=pd, color ="black", fill = "black", size=4)+
  scale_shape_manual(values=c(16,17),
                     name = "Origin") +
  annotate("text", x = 1, y = 0.25, label = "Rim > Patch") + 
  annotate("text", x = 1, y = -0.25, label = "Rim < Patch") + 
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank())+
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18, color = "black")) +
  theme(legend.position = c(0.9,0.8))

ggsave(file = "output/Graphs/L2018.vol.pdf", size.resid.2018)

# Statistics 
L2018.vol.comp$Origin.x<-droplevels(L2018.vol.comp$Origin.x)
L2018.vol.comp$Treatment.x<-droplevels(L2018.vol.comp$Treatment.x)

vol.resid.2018.anova <- lm(resid.1~Origin.x*Treatment.x, data = L2018.vol.comp)
qqnorm(resid(vol.resid.2018.anova))
qqline(resid(vol.resid.2018.anova)) 

boxplot(resid(vol.resid.2018.anova)~L2018.vol.comp$Origin.x)
boxplot(resid(vol.resid.2018.anova)~L2018.vol.comp$Treatment.x) #not normal

anova(vol.resid.2018.anova)

capture.output(anova(vol.resid.2018.anova), file = "output/Statistics/L2018.vol.csv")

