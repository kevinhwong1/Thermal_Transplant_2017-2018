#Title: Photosynthesis Irradiance Curves
#Author: HM Putnam, K Wong
#Date Last Modified: 20200723
#See Readme file for details

rm(list=ls()) #clears workspace 

#Read in required libraries
##### Include Versions of libraries

library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library("phytotools")

#PHOTOSYNTHESIS
path.p<-"data/2017/Colony_Respirometry/PI_Curve/" #the location of all your respirometry files
file.names<-list.files(path = path.p, pattern = "csv$") #list all csv file names in the folder
Photo.R <- data.frame(matrix(NA, nrow=length(file.names)*7, ncol=3)) #generate a 3 column dataframe with specific column names
colnames(Photo.R) <- c("Fragment.ID","Intercept", "µmol.L.sec")

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Photo.Data1 <-read.table(file.path(path.p,file.names[i]), skip = 1, header=T, sep=",", na.string="NA", fill = TRUE, as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Photo.Data1  <- Photo.Data1[,c(2,9,16)] #subset columns of interest
  Photo.0 <- subset(Photo.Data1, Time > "14:00:00" & Time < "14:20:00")
  Photo.15 <- subset(Photo.Data1, Time > "11:26:00" & Time < "11:59:00")
  Photo.30 <- subset(Photo.Data1, Time > "12:00:00" & Time < "12:18:00")
  Photo.50 <- subset(Photo.Data1, Time > "12:19:00" & Time < "12:41:00")
  Photo.65 <- subset(Photo.Data1, Time > "12:42:00" & Time < "13:03:00")
  Photo.85 <- subset(Photo.Data1, Time > "13:04:00" & Time < "13:22:00")
  Photo.100 <- subset(Photo.Data1, Time > "13:23:00" & Time < "13:40:00")
    
  lt.levs <- list(Photo.0,Photo.15,Photo.30,Photo.50,Photo.65,Photo.85,Photo.100)

  for(j in 1:length(lt.levs)){    
    Photo.Data <- as.data.frame(lt.levs[j])
  
    n<-dim(Photo.Data)[1] #identify length of data
    Photo.Data$sec <- 1:n #set seconds by one from start to finish of run
    
    #Save plot prior to and after data thinning to make sure thinning is not too extreme
    rename <- sub("_.*", "", file.names[i])
    pdf(paste0("output/Photo_Resp_Output/PI_curve/",rename,"_",j,"photo.thinning.pdf"))
    par(omi=rep(0.3, 4)) #set size of the outer margins in inches
    par(mfrow=c(1,2)) #set number of rows and columns in multi plot graphic
    plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot data as a function of time
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
    axis(1)
    axis(2, las=1)
    
    Photo.Data   <-  thinData(Photo.Data , by=20)$newData1 #thin data by every 20 points
    Photo.Data $sec <- as.numeric(rownames(Photo.Data )) #maintain numeric values for time
    plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot thinned data
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
    axis(1)
    axis(2, las=1)
    dev.off()
    
    Regs  <-  rankLocReg(xall=Photo.Data $sec, yall=Photo.Data $Value, alpha=0.2, 
                         method="pc", verbose=TRUE) 
    pdf(paste0("output/Photo_Resp_Output/PI_curve/",rename,"_",j,"photo.regs.pdf"))
    plot(Regs)
    dev.off()
    
    s <- seq(0,nrow(Photo.R),length(lt.levs)) #to order the file output sequence in correct order in data frame
    Photo.R[j+s[i],2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
    Photo.R[j+s[i],1] <- paste0(rename,"_",j) #stores the file name in the Date column
}
}

Photo.R

write.csv(Photo.R, file="data/2017/Colony_Respirometry/PI.data.csv")

#Load Sample Info
Sample.Info <- read.csv(file="data/2017/Colony_Respirometry/PI_Sample_Info.csv", header=T) #read in volume and sample.info data
Sample.Info$Vol.L <- Sample.Info$Chamber.Vol.mL/1000 #calculate volume
Sample.Info$Light <- Sample.Info$Light*1.32 # adjust light for Apogee correction

#Merge with sample info
Photo <- merge(Photo.R,Sample.Info, by="Fragment.ID")

#Load Surface Area standards from core area
Surf.Area <- 9.62 #cm2

#Account for blank rate
Blnk.Vol <- 0.620 #volume of a blank chamber in L

#Convert to µmol hr-1 
Photo$rate.hr <- (Photo$µmol.L.sec*60*60)*Photo$Vol.L #converts to µmol hr-1

#Subset by light level
L0 <- subset(Photo, Light.Level==0)
L15 <- subset(Photo, Light.Level==15)
L30 <- subset(Photo, Light.Level==30)
L50 <- subset(Photo, Light.Level==50)
L65 <- subset(Photo, Light.Level==65)
L85 <- subset(Photo, Light.Level==85)
L100 <- subset(Photo, Light.Level==100)

# Calculate blank chambers µmol hr-1
Amb.Blank.L0 <- (((aggregate(µmol.L.sec ~ Genotype, data=L0, mean))[1,2])*60*60)*Blnk.Vol
Amb.Blank.L15 <- (((aggregate(µmol.L.sec ~ Genotype, data=L15, mean))[1,2])*60*60)*Blnk.Vol
Amb.Blank.L30 <- (((aggregate(µmol.L.sec ~ Genotype, data=L30, mean))[1,2])*60*60)*Blnk.Vol
Amb.Blank.L50 <- (((aggregate(µmol.L.sec ~ Genotype, data=L50, mean))[1,2])*60*60)*Blnk.Vol
Amb.Blank.L65 <- (((aggregate(µmol.L.sec ~ Genotype, data=L65, mean))[1,2])*60*60)*Blnk.Vol
Amb.Blank.L85 <- (((aggregate(µmol.L.sec ~ Genotype, data=L85, mean))[1,2])*60*60)*Blnk.Vol
Amb.Blank.L100 <- (((aggregate(µmol.L.sec ~ Genotype, data=L100, mean))[1,2])*60*60)*Blnk.Vol

#Correct for Blank Chamber µmol hr-1
L0$corr.rate.hr <- L0$rate.hr-Amb.Blank.L0
L15$corr.rate.hr <- L15$rate.hr-Amb.Blank.L15
L30$corr.rate.hr <- L30$rate.hr-Amb.Blank.L30
L50$corr.rate.hr <- L50$rate.hr-Amb.Blank.L50
L65$corr.rate.hr <- L65$rate.hr-Amb.Blank.L65
L85$corr.rate.hr <- L85$rate.hr-Amb.Blank.L85
L100$corr.rate.hr <- L100$rate.hr-Amb.Blank.L100

#Convert to µmol cm-2 hr-1 (Normalize to Surface Area)
L0$Pnet_umol.cm2.hr <- L0$corr.rate.hr/Surf.Area
L15$Pnet_umol.cm2.hr <- L15$corr.rate.hr/Surf.Area
L30$Pnet_umol.cm2.hr <- L30$corr.rate.hr/Surf.Area
L50$Pnet_umol.cm2.hr <- L50$corr.rate.hr/Surf.Area
L65$Pnet_umol.cm2.hr <- L65$corr.rate.hr/Surf.Area
L85$Pnet_umol.cm2.hr <- L85$corr.rate.hr/Surf.Area
L100$Pnet_umol.cm2.hr <- L100$corr.rate.hr/Surf.Area

#remove blanks rows
L0 <- subset(L0, Genotype!= "Blank")
L15 <- subset(L15, Genotype!= "Blank")
L30 <- subset(L30, Genotype!= "Blank")
L50 <- subset(L50, Genotype!= "Blank")
L65 <- subset(L65, Genotype!= "Blank")
L85 <- subset(L85, Genotype!= "Blank")
L100 <- subset(L100, Genotype!= "Blank")

#combine data for all light levels into large data frame
PI <- rbind(L0, L15, L30, L50, L65, L85, L100)
PI$Fragment.ID <- sub("_.*", "", PI$Fragment.ID) #remove light steps 1-7 from fragment names

##### Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980) #####

#Plot curves
PAR <- as.numeric(PI$Light)
Pc <- as.numeric(PI$Pnet_umol.cm2.hr)

pdf("output/Graphs/PICurve.pdf")
par(mfrow=c(1,1))
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(-1, 2), cex.lab=0.8,cex.axis=0.8,cex=1, main="", adj = 0.05) #set plot info
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*cm^-2*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
curve.nlslrc.PA = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,
                      start=list(Am=(max(Pc)-min(Pc)),AQY=0.05,Rd=-min(Pc),theta=.001)) 

my.fit.PA <- summary(curve.nlslrc.PA ) #summary of model fit

#draw the curve using the model fit
curve.fitting.PA <- curve((1/(2*summary(curve.nlslrc.PA)$coef[4,1]))*(summary(curve.nlslrc.PA)$coef[2,1]*x+summary(curve.nlslrc.PA)$coef[1,1]-sqrt((summary(curve.nlslrc.PA)$coef[2,1]*x+summary(curve.nlslrc.PA)$coef[1,1])^2-4*summary(curve.nlslrc.PA)$coef[2,1]*summary(curve.nlslrc.PA)$coef[4,1]*summary(curve.nlslrc.PA)$coef[1,1]*x))-summary(curve.nlslrc.PA)$coef[3,1],lwd=2,col="blue",add=T)

#Extract the parameters

#Amax (max gross photosytnthetic rate) 
Pmax.gross <- my.fit.PA$parameters[1]

#AQY (quantum yield) alpha  
AQY <- my.fit.PA$parameters[2]

#Rd (dark respiration)
Rd <- my.fit.PA$parameters[3]

# Ik light saturation point
Ik <- Pmax.gross/AQY

# Ic light compensation point
Ic <- Rd/AQY

# Net photosynthetic rates
Pmax.net <- Pmax.gross - Rd

#output parameters into a table
PI.Output <- as.data.frame(rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic))
row.names(PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")

text(x=500, y=-0.25, labels="AQY =  0.015\nIk = 122.445\nIc = 56.892")

dev.off()






