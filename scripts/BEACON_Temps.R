#Title: Thermal Transplant 2017-2018 - BEACON field temps
#Author: KH Wong
#Date Last Modified: 202023
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

#Import data 
Crescent <- read.csv("data/2018/Field_Temp/Crescent_Reef_Aug2017_Jul2018.csv")
Hog <- read.csv("data/2018/Field_Temp/Hog_Reef_Aug2018_Jul2018.csv")

#Adding Reef Column
Crescent$Reef <- "Patch"
Hog$Reef <- "Rim"

#Binding datasets together
temp.all <- rbind(Crescent, Hog)

ggplot(temp.all, aes(x=Date, y = SST_C, group = Reef, color = Reef)) +
  geom_line()

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

#Temperature Data avg per day
Log_data_sum <- summarySE(Log_data, measurevar="Temp", groupvars=c("Date","Treatment")) #Summarizing by date and tank

ggplot(Log_data_sum,aes(x=Date, y=Temp, colour=Treatment))+
  geom_errorbar(aes(ymin=Temp-se, ymax=Temp+se), width=.1, colour="black") +
  ylab("Temperature")+
  scale_x_date(date_labels="%b %d",date_breaks  ="1 day") + #modifies how many dates are shown on the x axis
  geom_point()+
  geom_line(aes(group=Treatment))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Totals
ggplot(Log_data, aes(x=Tank, y=Temp))+ #boxplot of all data
  geom_boxplot() +
  ylab("Temperature")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Tank Stats
#Testing Assumptions
#Normality  
hist(Log_data$Temp)

logtemp<-log(Log_data$Temp)
hist(logtemp)

#Homogeneity of Variance
library(car)
leveneTest(Log_data$Temp~Log_data$Tank)
leveneTest(logtemp~Log_data$Tank)

#Nonparametric stats
kruskal.test(Temp~Tank, data=Log_data)
library(PMCMR)
posthoc.kruskal.nemenyi.test(Temp~Tank, data=Log_data,dist = "Tukey")


