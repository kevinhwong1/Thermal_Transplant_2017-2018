#Title: Thermal Transplant 2017-2018 - FECUNDITY
#Author: KH Wong
#Date Last Modified: 20181121
#See Readme file for details 

#Load Libaries
library(dbplyr)
library(tidyverse)
library(readr)
library(stringr)
library(lme4)
library(lmerTest)

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

fecundity <-  ggplot(fec.all, aes(x=treatment, y=total.plan.sa2 )) + #set up plot information
  geom_errorbar(aes(ymax=total.plan.sa2+se, ymin=total.plan.sa2-se), colour="black", width=.1, position = position_dodge(width = 0.2),stat = "identity") + #add standard error bars about the mean
  geom_point(aes(shape=reef.zone, color = Timepoint), position = position_dodge(width = 0.2), stat = "identity", size=4) + #plot points
#  scale_shape_manual(values=c(16,17)) + #sets point shape manually
#  scale_color_manual(values=c("blue", "red")) +
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


### Making Plasticity plots 

p

fec.tot.TPatch$coral.id2 <- fec.tot.TPatch$coral.id %>% str_replace("-A","")
fec.tot.TRim$coral.id2 <- fec.tot.TRim$coral.id %>% str_replace("-B","")

fec.tot.stat <- rbind(fec.tot.TPatch,fec.tot.TRim)
# drop extra blank factor
fec.tot.stat$reef.zone<-droplevels(fec.tot.stat$reef.zone)
fec.tot.stat$treatment<-droplevels(fec.tot.stat$treatment)

fec.2018.anova <- lmer(total.plan.sa2~reef.zone*treatment*Timepoint + (1|coral.id2), data = fec.tot.stat)
qqnorm(resid(fec.2018.anova))
qqline(resid(fec.2018.anova))

boxplot(resid(fec.2018.anova)~fec.tot.stat$reef.zone)
boxplot(resid(fec.2018.anova)~fec.tot.stat$treatment)
boxplot(resid(fec.2018.anova)~fec.tot.stat$Timepoint)
anova(fec.2018.anova)
summary(fec.2018.anova)


plasticity <- merge(fec.tot.TPatch, fec.tot.TRim, by = "coral.id2")
plasticity$resid.1 <- resid(lm(plasticity$total.plan.sa2.y - plasticity$total.plan.sa2.x ~0))


fec.resid.mean <- summarySE(plasticity, measurevar="resid.1", groupvars=c("treatment.x", "reef.zone.x"))

pd <- position_dodge(0.1) # moves object .05 to the left and right

resid.fecundity<- ggplot(fec.resid.mean, aes(x=treatment.x, y=resid.1, color=reef.zone.x, group=reef.zone.x)) + 
  geom_line(position=pd, color="black")+
  geom_errorbar(aes(ymin=resid.1-se, ymax=resid.1+se), width=.1, position=pd, color="black") + #Error bars
  geom_hline(yintercept = 0,linetype="dashed") +
#  ylim(0,16)+
  xlab("Treatment") + ylab(expression(paste("Mean Larval Release per" ~ cm^{2} , " (Rim - Patch)")))+ #Axis titles
  geom_point(position=pd, aes(fill=reef.zone.x), color ="black", pch=21, size=4)+
  scale_fill_discrete(name = "Origin") + 
  annotate("text", x = 1, y = 2, label = "Rim > Patch") + 
  annotate("text", x = 1, y = -10, label = "Rim < Patch") + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#all 
plot.plasticity <- ggplot(data = plasticity, aes(x=total.plan.sa2.x, y=total.plan.sa2.y))+
  ylab("Fecundity when transplanted to Rim")+ xlab("Fecundity when transplanted to Patch") + 
  ylim(0, 20)+ xlim(0,20) +
  geom_point(aes(shape=reef.zone.x, color = treatment.x), size = 3)+
  scale_shape_manual(values=c(16,17),
                     name = "Reef Zone")+
  scale_color_manual(values=c("blue", "red"),
                     name = "Treatment")+ #colour modification
  geom_abline(intercept = 0, slope = 1, linetype="dotted")+
  theme_bw() + theme(panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    panel.background = element_rect(colour = "black", size=1), legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Patch origin
patch.fec.plast <- plasticity %>%
  filter (reef.zone.x == "Patch")

plot.plasticity.patch <- ggplot(data = patch.fec.plast, aes(x=total.plan.sa2.x, y=total.plan.sa2.y))+
  ylab("Fecundity when transplanted to Rim")+ xlab("Fecundity when transplanted to Patch") + 
  ylim(0, 20)+ xlim(0,20) +
  geom_point(aes(shape=reef.zone.x, color = treatment.x), size = 3)+
  scale_shape_manual(values=c(16,17),
                     name = "Reef Zone")+
  scale_color_manual(values=c("blue", "red"),
                     name = "Treatment")+ #colour modification
  geom_abline(intercept = 0, slope = 1, linetype="dotted")+
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Rim origin
rim.fec.plast <- plasticity %>%
  filter (reef.zone.x == "Rim")

plot.plasticity.rim <- ggplot(data = rim.fec.plast, aes(x=total.plan.sa2.x, y=total.plan.sa2.y))+
  ylab("Fecundity when transplanted to Rim")+ xlab("Fecundity when transplanted to Patch") + 
  ylim(0, 20)+ xlim(0,20) +
  geom_point(aes(shape=reef.zone.x, color = treatment.x), size = 3)+
  scale_shape_manual(values=c(17, 16),
                     name = "Reef Zone")+
  scale_color_manual(values=c("blue", "red"),
                     name = "Treatment")+ #colour modification
  geom_abline(intercept = 0, slope = 1, linetype="dotted")+
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#fec plasticity graphs
# Larval Phys graph
fec.plastic.2018 <- arrangeGrob (plot.plasticity.patch, plot.plasticity.rim, ncol=2)
ggsave(file="~/MyProjects/Thermal_Transplant_2017-2018/output/fec.plastic.2018.pdf", fec.plastic.2018, width = 11, height = 6, units = c("in"))

