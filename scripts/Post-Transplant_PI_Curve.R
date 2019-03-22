#Title: PI curve for 2017-2018 Transplant Experiment 
#Author: HM Putnam and NJ Silbiger
#Edited by: Kevin Wong
#Date Last Modified: 20180719
#See Readme file for details 

##Install packages
if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 
library(devtools)
if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') install_github('colin-olito/LoLinR') 
if ("lubridate" %in% rownames(installed.packages()) == 'FALSE') install.packages('lubridate') 
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron') 
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 

#Read in required libraries

##### Include Versions of libraries
install_github('colin-olito/LoLinR')
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library('plyr')
library('dplyr')
library(stringr)

##PI CURVE FITTING

#PI.curve.data<-read.csv('~data/Colony.Respirometry/Transp.PI.Curve.csv', header=T, sep = ",")
install.packages('phytotools')
library('phytotools')

PAR <- c(0,
         33.6,
         59.2,
         106,
         152.2,
         287.4,
         473.6,
         602,
         803.6) #umol photons m-2 s-1
P.11.A <- c(-12.94590436,
            -9.045692,
            -7.924455866,
            -3.728947449,
            -0.962667146,
            1.035645869,
            3.503857592,
            4.019513747,
            2.576833866) #e.g., µmol O2 cm-2 sec-1

P.15.A <- c(-6.955172896,
            -3.762284871,
            -1.353798824,
            1.001380557,
            3.031387527,
            4.771539855,
            6.626304225,
            6.674334259,
            4.662077298) #e.g., µmol O2 cm-2 sec-1

R.4.A <- c(-7.053264844,
           -3.945091701,
           -1.18259728,
           1.941012164,
           4.303536245,
           7.137453234,
           9.688390597,
           9.276687659,
           10.89734604) #e.g., µmol O2 cm-2 sec-1
R.7.B <- c(-1.838393379,
           -2.703094952,
           -2.333213889,
           1.748332539,
           4.312864619,
           6.72017876,
           9.153658177,
           9.015297733) #e.g., µmol O2 cm-2 sec-1

#Model J.M. Heberling (2013) https://sites.google.com/site/fridleylab/home/protocols/sampleCode_lrc.R?attredirects=0
#https://academic.oup.com/aobpla/article/9/2/plx011/3074888#supplementary-data
curve.P.11.A = nls(P.11.A ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,start=list(Am=(max(P.11.A)-min(P.11.A)),AQY=0.05,Rd=-min(P.11.A),theta=1)) 

my.fit <- summary(curve.P.11.A) #summary of model fit
curve.fitting <- curve((1/(2*summary(curve.P.11.A)$coef[4,1]))*(summary(curve.P.11.A)$coef[2,1]*x+summary(curve.P.11.A)$coef[1,1]-sqrt((summary(curve.P.11.A)$coef[2,1]*x+summary(curve.P.11.A)$coef[1,1])^2-4*summary(curve.P.11.A)$coef[2,1]*summary(curve.P.11.A)$coef[4,1]*summary(curve.P.11.A)$coef[1,1]*x))-summary(curve.P.11.A)$coef[3,1],lwd=2,col="blue",add=T)

#Amax (max gross photosytnthetic rate) 
Pmax.gross <- my.fit$parameters[1]

#AQY (apparent quantum yield) alpha  
AQY <- my.fit$parameters[2]

#Rd (dark respiration)
Rd <- my.fit$parameters[3]

Ik <- Pmax.gross/AQY

Ic <- Rd/AQY

Pmax.net <- Pmax.gross - Rd

PI.Output <- rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic)
row.names(PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
PI.Output

#Plot input data and model fit
par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
plot(PAR,P.11.A,xlab="", ylab="", ylim=c(-15,max(P.11.A)+2),cex.lab=0.8,cex.axis=0.8,cex=1)
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1)
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*cm^-2*s^-1*")"),side=2,line=2,cex=1)
lines(curve.fitting ,lwd=2,col="blue")



