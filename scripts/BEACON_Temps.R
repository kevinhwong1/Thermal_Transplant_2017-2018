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

#Making a unique column
Crescent$Date.Time <- paste(Crescent$Date, Crescent$Time, sep = "-")
Hog$Date.Time <- paste(Hog$Date, Hog$Time, sep = "-")

#Merging data frames
meta.field.temp <- merge(Crescent, Hog, by = "Date.Time")

#Rename Columns



