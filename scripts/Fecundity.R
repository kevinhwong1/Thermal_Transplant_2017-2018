#Title: Thermal Transplant 2017-2018 - FECUNDITY
#Author: KH Wong
#Date Last Modified: 20181121
#See Readme file for details 

#Load Libaries
library(dbplyr)
library(tidyverse)
library(readr)
library(stringr)

#Set working directory
setwd("~/MyProjects/Thermal_Transplant_2017-2018/data/") 

### 2017 data ###

# load data 
pr.data2017<-read.csv('2017/TA_planrelease_2017.csv', header=T, sep=",")
coral.sa2017<-read.csv('2017/coral_sa.csv', header=T, sep = ",")

#Standardizing by SA for each vlaue 
pr.data.sa2017<-merge(pr.data2017, coral.sa2017, by="coral.id") #Merging data frames organized by coral.id
pr.data.sa2017 <- transform(pr.data.sa2017, total.plan.sa2=total.plan / sa) #divded total.plan by sa to standardize for sa and making a new variable

#Adding Timepoint Value
pr.data.sa2017$Timepoint <- "2"

#Summing per colony 
library(plyr)
fec.tot.2017<- ddply(pr.data.sa2017, .(coral.id, treatment,reef.zone,Timepoint), #summing rows by variables
                    summarize, 
                    total.plan.sa2 = sum(total.plan.sa2))

fec.tot.all.2 <- summarySE(fec.tot.2017, measurevar="total.plan.sa2", groupvars=c("treatment","reef.zone","Timepoint")) #Summarizing by treatment and reef


### 2018 data ###

pr.data2018<-read.csv('2018/Larval.Release.Transp.csv', header=TRUE, sep=",")
coral.sa2018<-read.csv('2018/coral.trans.calcsa.csv', header=TRUE, sep = ",")
info2018 <- read.csv('2018/Sample_Info_Transp.csv', header=TRUE, sep=",")

#renaming info columns
colnames(info2018)[colnames(info2018)=="Fragment.ID"] <- "coral.id"
colnames(info2018)[colnames(info2018)=="Treatment"] <- "treatment"
colnames(info2018)[colnames(info2018)=="Origin"] <- "reef.zone"

#Standardization and meta data attachment 
pr.data.sa2018 <- merge(pr.data2018, coral.sa2018, by="coral.id") #Merging data frames organizsed by coral.id
pr.data.sa2018 <- transform(pr.data.sa2018, total.plan.sa2=Number.Larvae / Surface_Area) #divded total.plan by sa to standardize for sa and making a new variable
pr.data.sa.2018 <- merge(pr.data.sa2018, info2018, by="coral.id")
pr.data.2018 <- pr.data.sa.2018 %>%
  select(coral.id, Date.x, Lunar.Day, total.plan.sa2, reef.zone, treatment, Transplant.Site)

#Adding timepoint values
pr.data.2018.TPatch <- pr.data.2018 %>%
  filter(Transplant.Site == "Patch")
pr.data.2018.TPatch$Timepoint <- "3P"

pr.data.2018.TRim <- pr.data.2018 %>%
  filter(Transplant.Site == "Rim")
pr.data.2018.TRim$Timepoint <- "3R"

# Summing data
fec.tot.TPatch<- ddply(pr.data.2018.TPatch, .(coral.id, treatment,reef.zone,Timepoint), #summing rows by variables
                     summarize, 
                     total.plan.sa2 = sum(total.plan.sa2))

fec.tot.all.TPatch <- summarySE(fec.tot.TPatch, measurevar="total.plan.sa2", groupvars=c("treatment","reef.zone","Timepoint")) #Summarizing by treatment and reef

fec.tot.TRim<- ddply(pr.data.2018.TRim, .(coral.id, treatment,reef.zone,Timepoint), #summing rows by variables
                       summarize, 
                       total.plan.sa2 = sum(total.plan.sa2))

fec.tot.all.TRim <- summarySE(fec.tot.TRim, measurevar="total.plan.sa2", groupvars=c("treatment","reef.zone","Timepoint")) #Summarizing by treatment and reef

# Merging ALL data sets

fec.all <- rbind(fec.tot.all.2, fec.tot.all.TPatch, fec.tot.all.TRim)
fec.all$Reef.treat <- paste(fec.all$reef.zone, fec.all$treatment)

# Plotting

fecundity <-  ggplot(fec.all, aes(x=Timepoint, y=total.plan.sa2, group=Reef.treat, color=treatment)) + #set up plot information
  geom_errorbar(aes(ymax=total.plan.sa2+se, ymin=total.plan.sa2-se), colour="black", width=.1, position = position_dodge(width = 0.2),stat = "identity") + #add standard error bars about the mean
  geom_point(aes(shape=reef.zone, color = treatment), position = position_dodge(width = 0.2), stat = "identity", size=4) + #plot points
  scale_shape_manual(values=c(16,17)) + #sets point shape manually
  scale_color_manual(values=c("blue", "red")) +
#  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  labs(x ="", y= "Mean larval release per cm2", treatment = "Treatment", reef.zone = "Origin") + #label y axis
  ylim(0, 17)+
#  scale_x_discrete(labels=c("2" = "Treatment Period (2017)", "3P" = "Transplanted to Patch (2018)",
#                            "3R" = "Transplanted to Rim (2018)")) +
  geom_vline(xintercept=1.5, linetype="dotted") +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_blank(), #Set the axes color
        panel.border = element_rect(linetype = "solid", color = "black"), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position=c(0.8, 0.75)) + #set legend location
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes
fecundity #view plot

ggsave(file="../Output/Fecundity.pdf", fecundity)
