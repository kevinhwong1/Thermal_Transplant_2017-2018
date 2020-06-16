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
library(plyr)
library(Rmisc)

### 2017 data ###

# load data 
pr.data2017<-read.csv('data/2017/Larval.Release/TA_planrelease_2017.csv', header=T, sep=",")
coral.sa2017<-read.csv('data/2017/Surface.Area/Colony.2017.sa.csv', header=T, sep = ",")

#Standardizing by SA for each vlaue 
pr.data.sa2017<-merge(pr.data2017, coral.sa2017, by="coral.id") #Merging data frames organized by coral.id
pr.data.sa2017 <- transform(pr.data.sa2017, total.plan.sa2=total.plan / sa) #divded total.plan by sa to standardize for sa and making a new variable

#Summing per colony 

fec.tot.2017<- ddply(pr.data.sa2017, .(coral.id, treatment,reef.zone), #summing rows by variables
                    summarize, 
                    total.plan.sa2 = sum(total.plan.sa2))

#Summarizing
fec.tot.2017.mean <- summarySE(fec.tot.2017, measurevar="total.plan.sa2", groupvars=c("treatment","reef.zone")) #Summarizing by treatment and reef

#Plotting
pd <- position_dodge(0.1) # moves object .05 to the left and right
legend.title <- "Treatment"
fec.2017<- ggplot(fec.tot.2017.mean, aes(x=treatment, y=total.plan.sa2, group=reef.zone)) + 
  geom_errorbar(aes(ymin=total.plan.sa2-se, ymax=total.plan.sa2+se), width=.1, position=pd, color="black") + #Error bars
  geom_line(position=pd, color="black", aes(linetype=reef.zone), size = 2)+
  ylim(-1,8)+
  xlab("Treatment") + ylab(expression("Mean Larval Release " ~ cm^{-2}))+ #Axis titles
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
fec.2017

ggsave(file = "output/Graphs/Fec.2017.point.pdf", fec.2017, width = 11, height = 11, units = c("in"))


#Violin plot
Violin.A2017.Fec <- ggplot(fec.tot.2017, aes(x=reef.zone, y=total.plan.sa2.log, fill = treatment)) +
  geom_violin(position = position_dodge(width = 0.9)) +
  geom_boxplot(width=.3, outlier.colour=NA, position = position_dodge(width = 0.9)) +
  geom_jitter(position = position_dodge(width = 0.9)) +
  #  stat_summary(fun.y=median, geom="line", position = position_dodge(width = 0.9), aes(group=Parental.Treatment))  + 
  #  stat_summary(fun.y=median, geom="point", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("#FFFFFF", "#999999")) +
  xlab("Reef Zone") + ylab(expression("Mean Larval Release " ~ cm^{-2})) + #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black")) +
  theme(legend.position = "none")

#ggsave(file = "output/Graphs/A2017.Zoox.violin.pdf", Violin.A2017.Zoox, width = 11, height = 11, units = c("in"))

p<-ggplot(fec.tot.2017, aes(x=total.plan.sa2, fill=treatment)) +
  geom_density(alpha=0.4) +
  scale_fill_manual(values=c("#FFFFFF", "#999999")) +
  facet_grid(reef.zone ~.)
  




# fec.2017.bar <- ggplot(fec.tot.2017.mean, aes(x=reef.zone, y=total.plan.sa2, fill=treatment)) + 
#   geom_bar(position=position_dodge(), stat="identity", color = "black") +
#   geom_errorbar(aes(ymin=total.plan.sa2-se, ymax=total.plan.sa2+se),
#                 width=.2,                    # Width of the error bars
#                 position=position_dodge(.9)) +
#   xlab("Reef Zone") + ylab(expression("Mean Larval Release per" ~ cm^{2}))+ #Axis titles
#   scale_fill_manual(values=c("dodgerblue3", "tomato1"),
#                     name = "Treatment") + #colour modification
#   theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
#                      panel.grid.minor = element_blank(), axis.line = element_blank()) +
#   theme(axis.text = element_text(size = 16, color = "black"),
#         axis.title = element_text(size = 18, color = "black"),
#         axis.title.x = element_blank()) +
#   theme(legend.position = c(0.9,0.8))
# 
# fec.2017.box <- ggplot(fec.tot.2017, aes(x=reef.zone, y=total.plan.sa2, fill=treatment)) + 
#   geom_boxplot()
# 
# 
#   xlab("Reef Zone") + ylab(expression("Mean Larval Release per" ~ cm^{2}))+ #Axis titles
#   scale_fill_manual(values=c("dodgerblue3", "tomato1"),
#                      name = "Treatment") + #colour modification
#   theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
#                      panel.grid.minor = element_blank(), axis.line = element_blank()) +
#   theme(axis.text = element_text(size = 16, color = "black"),
#         axis.title = element_text(size = 18, color = "black"),
#         axis.title.x = element_blank()) +
#   theme(legend.position = c(0.9,0.8))
# 
# ggsave(file = "output/Graphs/Fec.2017.bar.pdf", fec.2017.bar)

# Statistics model effects are robust to transformations --> sig to reef zone

fec.2017.anova1 <- lm(total.plan.sa2~reef.zone*treatment, data = fec.tot.2017)
qqnorm(resid(fec.2017.anova1)) #not normal
qqline(resid(fec.2017.anova1))
hist(resid(fec.2017.anova1))

boxplot(resid(fec.2017.anova)~fec.tot.2017$reef.zone)
boxplot(resid(fec.2017.anova)~fec.tot.2017$treatment)

anova(fec.2017.anova1)

#log
fec.tot.2017$total.plan.sa2.log <- log(fec.tot.2017$total.plan.sa2 + 1)

fec.2017.anova <- lm(total.plan.sa2.log~reef.zone*treatment, data = fec.tot.2017)
qqnorm(resid(fec.2017.anova)) #not normal
qqline(resid(fec.2017.anova))
hist(resid(fec.2017.anova))

boxplot(resid(fec.2017.anova)~fec.tot.2017$reef.zone)
boxplot(resid(fec.2017.anova)~fec.tot.2017$treatment)

anova(fec.2017.anova)

#sqrt
fec.tot.2017$total.plan.sa2.sqrt <- ((fec.tot.2017$total.plan.sa2 + 1)^(1/3))

fec.2017.anova.sqrt <- lm(total.plan.sa2.sqrt~reef.zone*treatment, data = fec.tot.2017)
qqnorm(resid(fec.2017.anova.sqrt)) #not normal
qqline(resid(fec.2017.anova.sqrt))
hist(resid(fec.2017.anova.sqrt))

boxplot(resid(fec.2017.anova)~fec.tot.2017$reef.zone)
boxplot(resid(fec.2017.anova)~fec.tot.2017$treatment)

anova(fec.2017.anova.sqrt)


#Removing variables with zero
fec.tot.2017.value <- fec.tot.2017 %>%
  filter(total.plan.sa2 > 0)

fec.tot.2017.value$total.plan.sa2.log <- log(fec.tot.2017.value$total.plan.sa2 + 1)

fec.2017.anova <- lm(total.plan.sa2.log~reef.zone*treatment, data = fec.tot.2017.value)
qqnorm(resid(fec.2017.anova)) #not normal
qqline(resid(fec.2017.anova))

boxplot(resid(fec.2017.anova)~fec.tot.2017.value$reef.zone)
boxplot(resid(fec.2017.anova)~fec.tot.2017.value$treatment)

anova(fec.2017.anova)

capture.output(anova(fec.2017.anova), file = "output/Statistics/A2017.Fec.csv")


######### Statistics 2020 analysis on 2017 fecundity data #########

## Question 1: What was the probability of releasing larvae based on reef site and treatment? 
# Binary logistic regression 





### 2018 data ###

pr.data2018<-read.csv('data/2018/Larval.Release/Larval.Release.Transp.csv', header=TRUE, sep=",")
coral.sa2018<-read.csv('data/2018/Surface.Area/colony.2018.calcsa.csv', header=TRUE, sep = ",")
info2018 <- read.csv('data/2018/Metadata/Sample_Info_Transp.csv', header=TRUE, sep=",")

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

fec.tot<- ddply(pr.data.2018, .(coral.id, treatment,reef.zone,Transplant.Site), #summing rows by variables
                       summarize, 
                       total.plan.sa2 = sum(total.plan.sa2))

fec.tot.mean <- summarySE(fec.tot, measurevar="total.plan.sa2", groupvars=c("treatment","reef.zone","Transplant.Site")) #Summarizing by treatment and reef
fec.tot.mean$reef.treatment <- paste(fec.tot.mean$reef.zone, fec.tot.mean$treatment)

pd <- position_dodge(0.1) # moves object .05 to the left and right
Fec2018 <- ggplot(fec.tot.mean, aes(x=Transplant.Site, y=total.plan.sa2, group=reef.treatment)) + 
  geom_line(position=pd, aes(linetype=reef.zone, color = treatment), size = 2)+
  geom_errorbar(aes(ymin=total.plan.sa2-se, ymax=total.plan.sa2+se), width=.1, position=pd, color="black") + #Error bars
  ylim(-2,16.1)+
  xlab("Transplant Site") + ylab(expression("Mean Larval Release " (cm^{-2})))+ #Axis titles
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

Fec2018
ggsave(file = "output/Graphs/A2018.Fec.pdf", Fec2018, width = 11, height = 11, units = c("in"))


# Statistics

fec.tot$reef.zone <- factor(fec.tot$reef.zone)
fec.tot$treatment <- factor(fec.tot$treatment)
fec.tot$Transplant.Site <- factor(fec.tot$Transplant.Site)

Fec2018Adult.anova <- lm(total.plan.sa2~reef.zone*treatment*Transplant.Site, data = fec.tot)
qqnorm(resid(Fec2018Adult.anova))
qqline(resid(Fec2018Adult.anova)) 

boxplot(resid(Fec2018Adult.anova)~fec.tot$reef.zone)
boxplot(resid(Fec2018Adult.anova)~fec.tot$treatment) 
boxplot(resid(Fec2018Adult.anova)~fec.tot$Transplant.Site)

anova(Fec2018Adult.anova)

capture.output(anova(Fec2018Adult.anova), file = "output/Statistics/A2018.Fec.csv")





# ### Residual Analysis ###
# #Adding timepoint values
# pr.data.2018.TPatch <- pr.data.2018 %>%
#   filter(Transplant.Site == "Patch")
# pr.data.2018.TPatch$Timepoint <- "3P"
# 
# pr.data.2018.TRim <- pr.data.2018 %>%
#   filter(Transplant.Site == "Rim")
# pr.data.2018.TRim$Timepoint <- "3R"
# 
# # Summing data
# fec.tot.TPatch<- ddply(pr.data.2018.TPatch, .(coral.id, treatment,reef.zone,Timepoint), #summing rows by variables
#                      summarize, 
#                      total.plan.sa2 = sum(total.plan.sa2))
# 
# fec.tot.all.TPatch <- summarySE(fec.tot.TPatch, measurevar="total.plan.sa2", groupvars=c("treatment","reef.zone","Timepoint")) #Summarizing by treatment and reef
# 
# fec.tot.TRim<- ddply(pr.data.2018.TRim, .(coral.id, treatment,reef.zone,Timepoint), #summing rows by variables
#                        summarize, 
#                        total.plan.sa2 = sum(total.plan.sa2))
# 
# fec.tot.all.TRim <- summarySE(fec.tot.TRim, measurevar="total.plan.sa2", groupvars=c("treatment","reef.zone","Timepoint")) #Summarizing by treatment and reef
# 
# #Making a new column with ust coral ID
# fec.tot.TPatch$coral.id2 <- fec.tot.TPatch$coral.id %>% str_replace("-A","")
# fec.tot.TRim$coral.id2 <- fec.tot.TRim$coral.id %>% str_replace("-B","")
# 
# #Re-binding Patch and Rim data frames
# fec.2018.comp <- merge(fec.tot.TPatch, fec.tot.TRim, by = "coral.id2")
# 
# #Making residuals between transplant sites
# fec.2018.comp$resid.1 <- resid(lm(fec.2018.comp$total.plan.sa2.y - fec.2018.comp$total.plan.sa2.x ~0))
# 
# #Summarizing residuals
# fec.resid.mean <- summarySE(fec.2018.comp, measurevar="resid.1", groupvars=c("treatment.x", "reef.zone.x"))
# 
# pd <- position_dodge(0.1) # moves object .05 to the left and right
# resid.fecundity<- ggplot(fec.resid.mean, aes(x=treatment.x, y=resid.1, group=reef.zone.x, shape = reef.zone.x)) + 
#   geom_line(position=pd, color="black", linetype = "3313")+
#   geom_errorbar(aes(ymin=resid.1-se, ymax=resid.1+se), width=.1, position=pd, color="black") + #Error bars
#   geom_hline(yintercept = 0,linetype="dashed") +
# #  ylim(0,16)+
#   xlab("Treatment") + ylab(expression(paste("Mean Larval Release per" ~ cm^{2} , " (Rim - Patch)")))+ #Axis titles
#   geom_point(aes(shape=reef.zone.x), position=pd, color ="black", fill = "black", size=4)+
#   scale_shape_manual(values=c(16,17),
#                      name = "Origin") +
#   annotate("text", x = 1, y = 2, label = "Rim > Patch") + 
#   annotate("text", x = 1, y = -10, label = "Rim < Patch") + 
#   theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
#                      panel.grid.minor = element_blank(), axis.line = element_blank())+
#   theme(axis.text = element_text(size = 16, color = "black"),
#         axis.title = element_text(size = 18, color = "black")) +
#   theme(legend.position = c(0.9,0.8))
# 
# ggsave(file = "output/Graphs/Fec.2018.resid.pdf", resid.fecundity)
# 
# 
# # drop extra blank factor
# fec.2018.comp$reef.zone.x<-droplevels(fec.2018.comp$reef.zone.x)
# fec.2018.comp$treatment.x<-droplevels(fec.2018.comp$treatment.x)
# 
# fec.resid.2018.anova <- lm(resid.1~reef.zone.x*treatment.x, data = fec.2018.comp)
# qqnorm(resid(fec.resid.2018.anova))
# qqline(resid(fec.resid.2018.anova)) #potentially not normal
# 
# boxplot(resid(fec.resid.2018.anova)~fec.2018.comp$reef.zone.x)
# boxplot(resid(fec.resid.2018.anova)~fec.2018.comp$treatment.x) #not normal
# 
# anova(fec.resid.2018.anova)
# summary(fec.resid.2018.anova)


