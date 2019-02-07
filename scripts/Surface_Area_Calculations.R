#Title: Thermal Transplant 2017-2018 - Coral Surface Area Calculations
#Author: KH Wong
#Date Last Modified: 20181121
#See Readme file for details 

library(ggplot2)
library(Rmisc)
library(dplyr)
library(tidyverse)

#Set working directory
setwd("~/MyProjects/Thermal_Transplant_2017-2018/data/") 

### 2018 coral surface area ###
coral.sa<-read.csv('2018/coral.trans.sa.csv', header=T, sep = ",")

#Subset standard
standard <- coral.sa %>% 
  filter(Type == "Standard") 

ggplot(data = standard, aes(x=Surface_Area, y=Weight))+
  ylab("Weight (g)")+ xlab("Surface Area (cm2)") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

lmstandard <- lm (Surface_Area ~ Weight, data = standard)
summary(lmstandard)

# Calculate SA with standard lm
coral.sa.calc <- coral.sa %>%
  filter(Type == "Coral")

coral.sa.calc$Surface_Area <- predict(lmstandard, newdata = coral.sa.calc) #using model to get concentration

write.csv(coral.sa.calc, '2018/coral.trans.calcsa.csv')
