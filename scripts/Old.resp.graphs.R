## PLasticity Plots

setwd("~/Dropbox/BIOS_Summer_2018/20172018_Transplant/RAnalysis/Output/") 


#Calculate means
AllMeans <- ddply(resp.data, c('Transplant.Site','Origin.Treatment', 'Origin', "Treatment"), summarize, #take out treatment?
                  #pnet
                  Pnet.mean= mean(Pnet_umol.cm2.hr, na.rm=T), #mean pnet
                  N = sum(!is.na(Pnet_umol.cm2.hr)), # sample size
                  Pnet.se = sd(Pnet_umol.cm2.hr, na.rm=T)/sqrt(N), #SE
                  #Rdark
                  Rdark.mean= mean(Rdark_umol.cm2.hr, na.rm=T), #mean rdark
                  Rdark.se = sd(Rdark_umol.cm2.hr, na.rm=T)/sqrt(N), #SE
                  #Pgross
                  Pgross.mean  = mean(Pgross_umol.cm2.hr, na.rm=TRUE),
                  Pgross.se = sd (Pgross_umol.cm2.hr, na.rm=TRUE)/sqrt(N))
AllMeans
#write Results
write.csv(resp.data, file="../Data/AllRates.csv") # raw data
write.csv(AllMeans, file="../Output/AllMeans.csv") # Mean data

#Plot the pnet mean
Fig.Pn <-  ggplot(AllMeans, aes(x=Transplant.Site, y=Pnet.mean,  group=Origin.Treatment)) + #set up plot information
  geom_errorbar(aes(x=Transplant.Site, ymax=Pnet.mean+Pnet.se, ymin=Pnet.mean-Pnet.se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=Origin, color = Treatment), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(16,17)) + #sets point shape manually
  scale_color_manual(values=c("blue", "red")) +
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Transplant Site") + #label x axis
  ylab("Net Photosynthesis umol cm-2 h-1") + #label y axis
  ylim(0, 1)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_blank(), #Set the axes color
        panel.border = element_rect(linetype = "solid", color = "black"), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position='none') + #set legend location
  ggtitle("Net Photosythesis") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes
Fig.Pn #view plot

#plot pgross mean
Fig.Pg <-  ggplot(AllMeans, aes(x=Transplant.Site, y=Pgross.mean, group=Origin.Treatment)) + #set up plot information
  geom_errorbar(aes(x=Transplant.Site, ymax=Pgross.mean+Pgross.se, ymin=Pgross.mean-Pgross.se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=Origin, color = Treatment), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(16,17)) + #sets point shape manually
  scale_color_manual(values=c("blue", "red")) +
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Transplant Site") + #label x axis
  ylab("Gross Photosynthesis umol cm-2 h-1") + #label y axis
  ylim(1,2.5)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_blank(), #Set the axes color
        panel.border = element_rect(linetype = "solid", color = "black"), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position='none') + #set legend location
  ggtitle("Gross Photosythesis") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes
Fig.Pg #view plot

#plot rdark mean
Fig.Rd <-  ggplot(AllMeans, aes(x=Transplant.Site, y=Rdark.mean,  group=Origin.Treatment)) + #set up plot information
  geom_errorbar(aes(x=Transplant.Site, ymax=Rdark.mean+Rdark.se, ymin=Rdark.mean-Rdark.se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=Origin, color = Treatment), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(16,17)) + #sets point shape manually
  scale_color_manual(values=c("blue", "red")) +
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Transplant Site") + #label x axis
  ylab("Dark Respiration umol cm-2 h-1") + #label y axis
  ylim(-1.25, 0)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_blank(), #Set the axes color
        panel.border = element_rect(linetype = "solid", color = "black"), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position=c(0.8, 0.75)) + #set legend location
  ggtitle("Dark Respiration") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes
Fig.Rd #view plot

dev.off()

Resp.Figs <- arrangeGrob( Fig.Pn, Fig.Pg,Fig.Rd, ncol=3)
ggsave(file="~/Dropbox/BIOS_Summer_2018/20172018_Transplant/RAnalysis/Output/Photo_Resp_Output/Transp_Respirometry.pdf", Resp.Figs, width = 11, height = 6, units = c("in"))

## Statsitics (time 2)
# Pnet
hist(resp.data2$Pnet_umol.cm2.hr)
library(car)
shapiro.test(resp.data2$Pnet_umol.cm2.hr) 
bartlett.test(Pnet_umol.cm2.hr~Treatment, data=resp.data2)
bartlett.test(Pnet_umol.cm2.hr~Origin, data=resp.data2)

Pnet.aov <- aov(Pnet_umol.cm2.hr ~ Treatment*Origin, data = resp.data2)
summary(Pnet.aov)
TukeyHSD(Pnet.aov)

#Pgross
hist(resp.data2$Pgross_umol.cm2.hr)
library(car)
shapiro.test(resp.data2$Pgross_umol.cm2.hr) 
bartlett.test(Pgross_umol.cm2.hr~Treatment, data=resp.data2)
bartlett.test(Pgross_umol.cm2.hr~Origin, data=resp.data2)

Pgross.aov <- aov(Pgross_umol.cm2.hr ~ Treatment*Origin, data = resp.data2)
summary(Pgross.aov)
TukeyHSD(Pgross.aov)

#Rdark

hist(resp.data2$Rdark_umol.cm2.hr)
library(car)
shapiro.test(resp.data2$Rdark_umol.cm2.hr) 
bartlett.test(Rdark_umol.cm2.hr~Treatment, data=resp.data2)
bartlett.test(Rdark_umol.cm2.hr~Origin, data=resp.data2)

Rdark.aov <- aov(Rdark_umol.cm2.hr ~ Treatment*Origin, data = resp.data2)
summary(Rdark.aov)
TukeyHSD(Rdark.aov)


## PLasticity Plots

setwd("~/Dropbox/BIOS_Summer_2018/20172018_Transplant/RAnalysis/Output/") 
resp.data <- read.csv("All.resp.data.csv")

resp.data.2018.TPatch <- resp.data %>%
  filter(Transplant.Site == "Patch")

resp.data.2018.TRim <- resp.data %>%
  filter(Transplant.Site == "Rim")

resp.data.2018.TPatch$coral.id <- resp.data.2018.TPatch$Fragment.ID %>% str_replace("-A","")
resp.data.2018.TRim$coral.id <- resp.data.2018.TRim$Fragment.ID %>% str_replace("-B","")

plasticity.resp <- merge(resp.data.2018.TPatch, resp.data.2018.TRim, by = "coral.id")

##residual plot Pnet
plasticity.resp$resid.1 <- resid(lm(plasticity.resp$Pnet_umol.cm2.hr.y - plasticity.resp$Pnet_umol.cm2.hr.x ~0))


resp.resid.mean <- summarySE(plasticity.resp, measurevar="resid.1", groupvars=c("Treatment.x", "Origin.x"))

pd <- position_dodge(0.1) # moves object .05 to the left and right
resid.resp<- ggplot(resp.resid.mean, aes(x=Treatment.x, y=resid.1, color=Origin.x, group=Origin.x)) + 
  geom_line(position=pd, color="black")+
  geom_errorbar(aes(ymin=resid.1-se, ymax=resid.1+se), width=.1, position=pd, color="black") + #Error bars
  geom_hline(yintercept = 0,linetype="dashed") +
  #  ylim(0,16)+
  xlab("Treatment") + ylab("Net Photosynthesis umol cm-2 h-1 (Rim - Patch)") + #Axis titles
  geom_point(position=pd, aes(fill=Origin.x), color ="black", pch=21, size=4)+
  scale_fill_discrete(name = "Origin") + 
  annotate("text", x = 1, y = 1, label = "Rim > Patch") + 
  annotate("text", x = 1, y = -0.5, label = "Rim < Patch") + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

Pnet.2018.anova <- lm(resid.1~Origin.x*Treatment.x, data = plasticity.resp)
qqnorm(resid(Pnet.2018.anova))
qqline(resid(Pnet.2018.anova))

boxplot(resid(Pnet.2018.anova)~plasticity.resp$Origin.x)
boxplot(resid(Pnet.2018.anova)~plasticity.resp$Treatment.x)

anova(Pnet.2018.anova)
summary(Pnet.2018.anova)

### residual plot resp
plasticity.resp$resid.2 <- resid(lm(plasticity.resp$Rdark_umol.cm2.hr.y - plasticity.resp$Rdark_umol.cm2.hr.x ~0))


resp.resid.mean2 <- summarySE(plasticity.resp, measurevar="resid.2", groupvars=c("Treatment.x", "Origin.x"))

pd <- position_dodge(0.1) # moves object .05 to the left and right
resid.resp2<- ggplot(resp.resid.mean, aes(x=Treatment.x, y=resid.2, color=Origin.x, group=Origin.x)) + 
  geom_line(position=pd, color="black")+
  geom_errorbar(aes(ymin=resid.2-se, ymax=resid.2+se), width=.1, position=pd, color="black") + #Error bars
  geom_hline(yintercept = 0,linetype="dashed") +
  #  ylim(0,16)+
  xlab("Treatment") + ylab("Dark Respiration umol cm-2 h-1 (Rim - Patch)") + #Axis titles
  geom_point(position=pd, aes(fill=Origin.x), color ="black", pch=21, size=4)+
  scale_fill_discrete(name = "Origin") + 
  annotate("text", x = 1, y = 0.5, label = "Rim > Patch") + 
  annotate("text", x = 1, y = -0.5, label = "Rim < Patch") + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



#all 
resp.plasticity <- ggplot(data = plasticity.resp, aes(x=Pnet_umol.cm2.hr.x, y=Pnet_umol.cm2.hr.y))+
  ylab("Oxygen Flux when transplanted to Rim")+ xlab("Oxygen Flux when transplanted to Patch") + 
  ylim(0, 1.2)+ xlim(0,1.2) +
  geom_point(aes(shape=Origin.x, color = Treatment.x), size = 3)+
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
patch.plasticity.resp <- plasticity.resp %>%
  filter (Origin.x == "Patch")

patch.resp.plasticity <- ggplot(data = patch.plasticity.resp, aes(x=Pnet_umol.cm2.hr.x, y=Pnet_umol.cm2.hr.y))+
  ylab("Oxygen Flux when transplanted to Rim")+ xlab("Oxygen Flux when transplanted to Patch") + 
  ylim(0, 1.2)+ xlim(0,1.2) +
  geom_point(aes(shape=Origin.x, color = Treatment.x), size = 3)+
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
rim.plasticity.resp  <- plasticity.resp %>%
  filter (Origin.x == "Rim")

rim.resp.plasticity <- ggplot(data = rim.plasticity.resp, aes(x=Pnet_umol.cm2.hr.x, y=Pnet_umol.cm2.hr.y))+
  ylab("Oxygen Flux when transplanted to Rim")+ xlab("Oxygen Flux when transplanted to Patch") + 
  ylim(0, 1.2)+ xlim(0,1.2) +
  geom_point(aes(shape=Origin.x, color = Treatment.x), size = 3)+
  scale_shape_manual(values=c(17, 16),
                     name = "Reef Zone")+
  scale_color_manual(values=c("blue", "red"),
                     name = "Treatment")+ #colour modification
  geom_abline(intercept = 0, slope = 1, linetype="dotted")+
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Resp plasticity graphs

resp.plastic.2018 <- arrangeGrob (patch.resp.plasticity, rim.resp.plasticity, ncol=2)
ggsave(file="~/MyProjects/Thermal_Transplant_2017-2018/output/resp.plastic.2018.pdf", resp.plastic.2018, width = 11, height = 6, units = c("in"))


