#Title: Thermal Transplant 2017-2018 - Coral Surface Area Calculations
#Author: KH Wong
#Date Last Modified: 20190320
#See Readme file for details 

library(ggplot2)
library(Rmisc)
library(dplyr)
library(tidyverse)


### 2017 colony surface area ###

### 2018 colony surface area ###
coral.sa.2018 <-read.csv("data/2018/Surface.Area/coral.trans.sa.csv", header=T, sep = ",")

#Subset standard
standard.2018 <- coral.sa.2018 %>% 
  filter(Type == "Standard") 

ggplot(data = standard.2018, aes(x=Surface_Area, y=Weight))+
  ylab("Weight (g)")+ xlab("Surface Area (cm2)") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

lmstandard.2018 <- lm (Surface_Area ~ Weight, data = standard.2018)
summary(lmstandard.2018)

# Calculate SA with standard lm
coral.sa.calc.2018 <- coral.sa.2018 %>%
  filter(Type == "Coral")

coral.sa.calc.2018$Surface_Area <- predict(lmstandard.2018, newdata = coral.sa.calc.2018) #using model to get concentration

write.csv(coral.sa.calc.2018, 'data/2018/Surface.Area/colony.2018.calcsa.csv')


### 2017 Fragment Surface Area ###
# Standard curve

frag.2017.standard <- read.csv("data/2017/Surface.Area/frag.SA.curve.csv") #loading data

plot.frag.2017.standard<- ggplot(data = frag.2017.standard, aes(x=Surface.Area.cm2, y=Mass.g))+
  ylab("Mass (g)")+ xlab("Surface Area (cm2)") + 
  geom_point()+
  geom_smooth(method = "lm") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


lm.frag.2017.standard <- lm (Surface.Area.cm2 ~ Mass.g, data = frag.2017.standard) #creating a linear model off of standard curve
lm.frag.2017.summary <- summary(lm.frag.2017.standard) #creating a summary outfut for the linear model
lm.frag.2017.summary 

# Calculating fragment surface area from weight

coral.frag.2017 <- read.csv("data/2017/Surface.Area/Adult.Frag.Processing.csv")
coral.frag.2017$Surface.Area <- predict(lm.frag.2017.standard, newdata = coral.frag.2017) #using model to get surface area
write.csv(coral.frag.2017, 'data/2017/Surface.Area/Adult.Frag.2017.Calculated.csv')

### 2018 Fragment Surface Area ### 


