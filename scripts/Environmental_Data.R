#Title: Thermal Transplant 2017-2018 - BEACON field temps
#Author: KH Wong
#Date Last Modified: 20200623
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
library(doBy)
library(devtools)
devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)
library(ggradar)
library(ggiraphExtra)
library(tidyverse)
library(scales)
library(tidyr)

##### 2017-2018 field temps #####
#Import data 
Crescent <- read.csv("data/2018/Field_Temp/Crescent_Reef_Aug2017_Jul2018.csv")
Hog <- read.csv("data/2018/Field_Temp/Hog_Reef_Aug2018_Jul2018.csv")

#Adding Reef Column
Crescent$Reef <- "Patch"
Hog$Reef <- "Rim"

#Binding datasets together
temp.all <- rbind(Crescent, Hog)
temp.all$Date <- as.Date(temp.all$Date, "%m/%d/%Y")
temp.all <- na.omit(temp.all)
temp.subset <- temp.all %>% # removing extra rim datapoints 
  filter(Date < "2018-06-10") %>%
  filter(Date > "2017-08-27")

raw.2017.2018 <- ggplot(temp.subset, aes(x=Date, y = SST_C, group = Reef, color = Reef)) +
  geom_point() + 
                    geom_line() +
  ylab("Temperature °C")+
                    scale_color_manual(values=c("tomato3", "dodgerblue3")) +
                    scale_x_date(date_breaks = "2 month", date_labels= "%b %Y") +
  scale_y_continuous(limits = c(17, 32), breaks = seq(17, 32, by = 2)) +
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1)) +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/2017.2018.temp.pdf", raw.2017.2018, width = 20, height = 11, units = c("in"))

# Summary statistics 
library(doBy)
temp.summary<- summaryBy(SST_C ~ Reef, data = temp.subset,
                         FUN = function(x) { c(max = max(x), min = min(x)) } )
temp.summary$Temp.seasonal <- temp.summary$SST_C.max - temp.summary$SST_C.min

# Daily Range
Daily.range <- summaryBy(SST_C ~ Date + Reef, data=temp.subset, FUN=c(max,min,mean,sd))

Daily.range$temp.daily.range <- Daily.range$SST_C.max - Daily.range$SST_C.min

temp.d.range <- summaryBy(temp.daily.range ~ Reef, data = Daily.range,
                          FUN = function(x) { c(mean = mean(x)) } )

temp.d.30 <- summaryBy(SST_C.max ~ Site, data = Daily.range,
                       FUN = function(x) { c(d.above.30 = count(x>30)) } )




##### 2011-2013 Temp and light #####

#Import data (Hog and Cresecent)
BDA.2011.2013 <- read.csv("data/2018/Field_Temp/BDA_2011-2013_TempLight.csv") 

#Changing date format 
BDA.2011.2013$Date <- as.Date(BDA.2011.2013$Date, "%m/%d/%Y")

#Removing rows with -999 values 
BDA.2011.2013.clean <- BDA.2011.2013 %>% 
  filter(Temp != -9999) %>%
  filter(Light != -9999) %>%
  filter(Temp < 31)  %>%
  filter(Date < "2012-09-05") %>%
  filter(Date > " 2010-10-18")

#Import data (Whalebone)

BDA.WB <- read.csv("data/2018/Field_Temp/Whalebone_2010-2012_Temp.csv") 

#Separating date and time column 
BDA.WB.2 <- BDA.WB %>%
  separate(Date_Time, c('Date', 'Time'), " ")

#Changing date format: Temp loggers switch format between 2011 and 2012, therefore we have to convert the dates separately then recombine
BDA.WB.2.a <- BDA.WB.2[1:35040,]
BDA.WB.2.b <- BDA.WB.2[35041:52608,]
BDA.WB.2.c <- BDA.WB.2[52609:87648,]
BDA.WB.2.d <- BDA.WB.2[87648:105216,]

BDA.WB.2.a$Date <- as.Date(BDA.WB.2.a$Date, "%m-%d-%y")
BDA.WB.2.c$Date <- as.Date(BDA.WB.2.c$Date, "%m-%d-%y")

BDA.WB.2.b$Date <- as.Date(BDA.WB.2.b$Date, "%Y-%m-%d")
BDA.WB.2.d$Date <- as.Date(BDA.WB.2.d$Date, "%Y-%m-%d")

BDA.WB.3 <- rbind(BDA.WB.2.a, BDA.WB.2.b, BDA.WB.2.c, BDA.WB.2.d)

# Removing NAs
BDA.WB.4 <- BDA.WB.3 %>%
  filter(Temp != "NA") %>%
  filter(Date < "2012-09-05") %>%
  filter(Date > " 2010-10-18")

#combining datasets
BDA.env.2010.2012.raw <- rbind(BDA.2011.2013.clean[,c(1,4,5,6,7)], BDA.WB.4)

#Averaging loggers
BDA.env.2010.2012 <- summarySE(BDA.env.2010.2012.raw, measurevar="Temp", groupvars=c("Site","Date","Time")) 

write.csv(BDA.env.2010.2012, file = 'data/2018/Field_Temp/All.Temp.data.csv')
# 
# # Mean Temp by date
# mean.temps.2011.2013 <- summarySE(BDA.2011.2013.clean, measurevar="Temp", groupvars=c("Date","Site"))
# mean.temps.2011.2013.select <- mean.temps.2011.2013 %>%
#   filter(Date < "2012-09-05") %>%
#   filter(Date > " 2010-10-18")
# 
# Mean.temp.20112013.plot <- ggplot(mean.temps.2011.2013.select, aes(x=Date, y = Temp, group = Site, color = Site)) + 
#                             geom_point() + 
#                             geom_line() +
# #                            geom_ribbon(aes(ymin=(Temp - ci), ymax=(Temp + ci)), linetype=2, alpha=0.1) +
#                             ylab("Temperature °C")+
#                             scale_x_date(date_breaks = "3 month", date_labels= "%b %Y") +
#                             scale_y_continuous(limits = c(17, 32), breaks = seq(17, 32, by = 2)) +
#                             scale_shape_manual(values=c(21,24),
#                                                name = "Site")+
#                             scale_color_manual(values=c("tomato3", "dodgerblue3")) +
#   theme_bw() + theme(panel.grid.major = element_blank(), 
#                      panel.grid.minor = element_blank(),
#                      panel.background = element_rect(colour = "black", size=1)) +
#   theme(axis.text = element_text(size = 30, color = "black"),
#         axis.title = element_text(size = 36, color = "black"),
#         axis.title.x = element_blank()) +
#   theme(legend.position = "none")
# 
# ggsave(file = "output/Graphs/2010.2012.temp.pdf", Mean.temp.20112013.plot, width = 20, height = 11, units = c("in"))

# Summary statistics 
library(doBy)
temp.summary<- summaryBy(Temp ~ Site, data = BDA.env.2010.2012,
          FUN = function(x) { c(max = max(x), min = min(x)) } )
temp.summary$Temp.seasonal <- temp.summary$Temp.max - temp.summary$Temp.min

# Daily Range
Daily.range <- summaryBy(Temp ~ Date + Site, data=BDA.env.2010.2012, FUN=c(max,min,mean,sd))

Daily.range$temp.daily.range <- Daily.range$Temp.max - Daily.range$Temp.min

Daily.range.clean <- Daily.range %>%
  filter(temp.daily.range < 4)

temp.d.range <- summaryBy(temp.daily.range ~ Site, data = Daily.range.clean,
                          FUN = function(x) { c(mean = mean(x)) } )

temp.d.30 <- summaryBy(Temp.max ~ Site, data = Daily.range.clean,
                          FUN = function(x) { c(d.above.30 = count(x>30)) } )
Temp.days.30 <- c(48, 4, 6)

temp.summary2 <- merge(temp.summary, temp.d.range, by = "Site")
temp.summary3 <- data.frame(temp.summary2, Temp.days.30)

temp.summary.4 <- temp.summary3 %>%
  filter(Site != "Whalebone")

temp.summary.5 <- temp.summary.4[2:6]

# light.summary<- summaryBy(Light ~ Site, data = BDA.2011.2013.clean,
#                           FUN = function(x) { c(m = mean(x), var = sd(x), max = max(x)) } )
# 
# temp.light.summary <- merge(temp.summary3, light.summary, by = "Site")
# temp.light.summary2 <- temp.light.summary[c(2:9)]


#### RADAR PLOT ####

fn <- function(x) x * 100/max(x, na.rm = TRUE)
fn(c(0,1,0))

## to all columns of your data frame
rad.data <- data.frame(lapply(temp.summary.5, fn))
rownames(rad.data) <- c("Patch", "Rim")
colnames(rad.data) <- c("Max", "Min", "Seasonal", "Daily", "30°C")
rad.data <- rad.data[,c(1,2,3,4,5)]

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
rad.data <- rbind(rep(100,5) , rep(0,5) , rad.data)

# Color vector
colors_border=c("tomato3", "dodgerblue3")
colors_in=c("tomato3", "dodgerblue3")

# Color vector
colors_border=c( rgb(0.9,0,0,0.9), rgb(0,0,0.9,0.9) )
colors_in=c( rgb(0.9,0,0,0.4), rgb(0,0,0.9,0.4))

library(fmsb)
# plot with default options:
radarchart(rad.data, axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,100,25), cglwd=0.8,
            #custom labels
            vlcex=0.9
)

dev.off()

# Add a legend
legend(x=0.9, y=0.9, legend = rownames(rad.data[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "black", cex=1.2, pt.cex=3)





###### EXPERIMENTAL TEMPERATURE 2017 ########


#Tank 1: Patch Ambient, Tank 2: Rim Ambient, Tank 3: Patch Hot, Tank 4: Rim HoT

### Data Manipulation

#load data (Acclimation)
PA <-read.csv('data/2017/Tank.Measurements/Logger/PA-705.csv', header=T, sep=",")
PB <-read.csv('data/2017/Tank.Measurements/Logger/RB-705.csv', header=T, sep=",")
RA <-read.csv('data/2017/Tank.Measurements/Logger/RA-705.csv', header=T, sep=",")
RB <-read.csv('data/2017/Tank.Measurements/Logger/RB-705.csv', header=T, sep=",")

PA$Tank <- rep(1,nrow(PA)) #adding a row with corresponding tank number
PB$Tank <- rep(1,nrow(PB))
RA$Tank <- rep(2,nrow(RA))
RB$Tank <- rep(2,nrow(RB))

PA$Treatment <- rep("Ambient",nrow(PA)) #adding a column with corresponding tank treatment
PB$Treatment <- rep("Ambient",nrow(PB))
RA$Treatment <- rep("Ambient",nrow(RA))
RB$Treatment <- rep("Ambient",nrow(RB))

# load data (Treatment JUL 5-11)
T1A <-read.csv('data/2017/Tank.Measurements/Logger/T1A-711.csv', header=T, sep=",")
T1B <-read.csv('data/2017/Tank.Measurements/Logger/T1B-711.csv', header=T, sep=",")
T2A <-read.csv('data/2017/Tank.Measurements/Logger/T2A-711.csv', header=T, sep=",")
T2B <-read.csv('data/2017/Tank.Measurements/Logger/T2B-711.csv', header=T, sep=",")
T3A <-read.csv('data/2017/Tank.Measurements/Logger/T3A-711.csv', header=T, sep=",")
T3B <-read.csv('data/2017/Tank.Measurements/Logger/T3B-711.csv', header=T, sep=",")
T4A <-read.csv('data/2017/Tank.Measurements/Logger/T4A-711.csv', header=T, sep=",")
T4B <-read.csv('data/2017/Tank.Measurements/Logger/T4B-711.csv', header=T, sep=",")

T1A$Tank <- rep(1,nrow(T1A)) #adding a coumn with corresponding tank number
T1B$Tank <- rep(1,nrow(T1B))
T2A$Tank <- rep(2,nrow(T2A))
T2B$Tank <- rep(2,nrow(T2B))
T3A$Tank <- rep(3,nrow(T3A))
T3B$Tank <- rep(3,nrow(T3B))
T4A$Tank <- rep(4,nrow(T4A))
T4B$Tank <- rep(4,nrow(T4B))

T1A$Treatment <- rep("Ambient",nrow(T1A)) #adding a column with corresponding tank treatment
T1B$Treatment <- rep("Ambient",nrow(T1B))
T2A$Treatment <- rep("Ambient",nrow(T2A))
T2B$Treatment <- rep("Ambient",nrow(T2B))
T3A$Treatment <- rep("Heated",nrow(T3A))
T3B$Treatment <- rep("Heated",nrow(T3B))
T4A$Treatment <- rep("Heated",nrow(T4A))
T4B$Treatment <- rep("Heated",nrow(T4B))

# load data (Treatment JUL 12-30)
Log_data_Treat<-read.csv('data/2017/Tank.Measurements/Logger/Logger_Data_Sum_2.csv', header=T, sep=",")

T1C <- Log_data_Treat %>% 
  filter(Tank == "1")

T2C <- Log_data_Treat %>% 
  filter(Tank == "2")

T3C <- Log_data_Treat %>% 
  filter(Tank == "3")

T4C <- Log_data_Treat %>% 
  filter(Tank == "4")

T1C$Treatment <- rep("Ambient",nrow(T1C)) #adding a column with corresponding tank treatment
T2C$Treatment <- rep("Ambient",nrow(T2C))
T3C$Treatment <- rep("Heated",nrow(T3C))
T4C$Treatment <- rep("Heated",nrow(T4C))

Log_data <- rbind(PA, PB, RA, RB, T1A, T1B, T2A, T2B, T3A, T3B, T4A, T4B, T1C, T2C, T3C, T4C) #Binding all data sets together

#date modification
Log_data$Date <- as.Date(Log_data$Date, format= "%d-%b-%Y")

Log_data$date.time<- paste(Log_data$Date, Log_data$Time, sep=" ") #combining date and time columns together

Log_data$Tank<- as.factor(Log_data$Tank) #Converting tank number from a numeric to a variable factor

## Subsetting acclimation and treatment 

acclim <- Log_data %>%
  filter(Date < "0017-07-05")

treatment <- Log_data %>%
  filter(Date > "0017-07-05")

#Temperature Data avg per day
Log_data_sum <- summarySE(Log_data, measurevar="Temp", groupvars=c("Date","Treatment")) #Summarizing by date and tank

tank.treat.time <- ggplot(Log_data_sum,aes(x=Date, y=Temp, colour=Treatment))+
#  geom_errorbar(aes(ymin=Temp-se, ymax=Temp+se), width=.1, colour="black") +
  ylab("Temperature °C")+
  scale_x_date(date_labels="%b %d",date_breaks  ="5 day") + #modifies how many dates are shown on the x axis
  scale_y_continuous(limits = c(28, 32), breaks = seq(28, 32, by = 1)) +
#  geom_point(size = 12)+
  geom_line(aes(group=Treatment), size = 3)+
  geom_ribbon(aes(ymin=(Temp - ci), ymax=(Temp + ci)), linetype=2, alpha=0.1) +
  scale_fill_manual(values=c("#FFFFFF", "#999999"),
                    name = "Treatment")+ #colour modification
  scale_color_manual(values=c("black", "#999999"),
                     name = "Treatment")+
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1)) +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")


ggsave(file = "output/Graphs/Treatment.Time.temp.pdf", tank.treat.time, width = 15, height = 11, units = c("in"))

#Totals
treatment.box <- ggplot(treatment, aes(x=Treatment, y=Temp, fill = Treatment))+ #boxplot of all data
  geom_boxplot() +
  ylab("Temperature °C")+
  scale_y_continuous(limits = c(28, 32), breaks = seq(27, 34, by = 1)) +
  scale_fill_manual(values=c("#FFFFFF", "#999999"),
                    name = "Treatment")+ #colour modification
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1)) +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/Treatment.box.temp.pdf", treatment.box, width = 11, height = 11, units = c("in"))



#Tank Stats
#Testing Assumptions
###ACCLIMATION###
#Normality  
hist(acclim$Temp)

logtemp<-log(acclim$Temp)
hist(logtemp)

#Homogeneity of Variance
library(car)
leveneTest(treatment$Temp~treatment$Tank)
leveneTest(logtemp~Log_data$Tank)

# Wilcox.test 
wilcox.test(Temp~Tank, data=acclim)

###TREATMENT###
#Normality  
hist(treatment$Temp)

logtemp<-log(treatment$Temp)
hist(logtemp)

#Homogeneity of Variance
library(car)
leveneTest(treatment$Temp~treatment$Tank)
leveneTest(logtemp~Log_data$Tank)

#Nonparametric stats
kruskal.test(Temp~Tank, data=treatment)
library(PMCMR)
posthoc.kruskal.nemenyi.test(Temp~Tank, data=treatment,dist = "Tukey")

# Wilcox.test 
wilcox.test(Temp~Treatment, data=treatment)


#Summary statistics 
Treatment.summary <- summarySE(treatment, measurevar="Temp", groupvars="Treatment")

Acc.summary <- summarySE(acclim, measurevar="Temp", groupvars="Treatment")

### Saility, pH, PAR

tank.env <- read.csv("data/2017/Tank.Measurements/Daily_Tank_Measures_Outside.csv")

tank.env.acc <- tank.env %>%
  filter(Experiment.Stage == "Acclimation")

acc.env.sal <- summarySE(tank.env.acc, measurevar="Salinity", groupvars="Experiment.Stage")
acc.env.pH <- summarySE(tank.env.acc, measurevar="pH", groupvars="Experiment.Stage")


tank.env.treat <- tank.env %>%
  filter(Experiment.Stage == "Treatment")

tank.1 <- tank.env.treat %>%
  filter(Tank == "Tank1")

tank.2 <- tank.env.treat %>%
  filter(Tank == "Tank2")

tank.3 <- tank.env.treat %>%
  filter(Tank == "Tank3")

tank.4 <- tank.env.treat %>%
  filter(Tank == "Tank4")

tank.1$Treatment <- "Ambient"
tank.2$Treatment <- "Ambient"
tank.3$Treatment <- "Hot"
tank.4$Treatment <- "Hot"

tank.env.treat.all <- rbind(tank.1, tank.2, tank.3, tank.4)
tank.env.treat.all$Treatment <- as.factor(tank.env.treat.all$Treatment)
tank.env.treat.all$Light.PAR.corr <- tank.env.treat.all$Light.PAR * 1.32

treat.env.sal <- summarySE(tank.env.treat.all, measurevar="Salinity", groupvars="Treatment")
treat.env.pH <- summarySE(tank.env.treat.all, measurevar="pH", groupvars="Treatment")
treat.env.PAR <- summarySE(tank.env.treat.all, measurevar="Light.PAR.corr", groupvars="Treatment")


#Stats for Salinity
sal.lm <- lm(Salinity~Treatment, data = tank.env.treat.all)
qqnorm(resid(sal.lm))
qqline(resid(sal.lm))

boxplot(resid(sal.lm)~tank.env.treat.all$Treatment)

t.test(Salinity~Treatment, data = tank.env.treat.all)

capture.output(t.test(Salinity~Treatment, data = tank.env.treat.all), file = "output/Statistics/salinity.csv")

#Stats for pH
pH.lm <- lm(pH~Treatment, data = tank.env.treat.all)
qqnorm(resid(pH.lm))
qqline(resid(pH.lm))

boxplot(resid(pH.lm)~tank.env.treat.all$Treatment)

t.test(pH~Treatment, data = tank.env.treat.all)

capture.output(t.test(Salinity~Treatment, data = tank.env.treat.all), file = "output/Statistics/salinity.csv")

#Stats for light
light.lm <- lm(Light.PAR.corr~Treatment, data = tank.env.treat.all)
qqnorm(resid(light.lm))
qqline(resid(light.lm))

boxplot(resid(light.lm)~tank.env.treat.all$Treatment)

t.test(Light.PAR.corr~Treatment, data = tank.env.treat.all)

capture.output(t.test(Salinity~Treatment, data = tank.env.treat.all), file = "output/Statistics/salinity.csv")

