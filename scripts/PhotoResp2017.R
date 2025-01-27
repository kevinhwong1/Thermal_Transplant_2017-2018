#Title: Photosynthesis and Respiration Calculations for Thermal Acclim. Proj 2017
#Author: HM Putnam and NJ Silbiger
#Edited by: Kevin Wong
#Date Last Modified: 20200716
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
library('Rmisc')


##### PHOTOSYNTHESIS Time 0 #####
path.p<-"data/2017/Colony_Respirometry/T1_clay/" #the location of all your respirometry files

# bring in the respiration files
file.names<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders
#basename above removes the subdirectory name from the file

#add file names that include the subdirectory name
file.names.full<-list.files(path = path.p, pattern = "csv$", recursive = TRUE) #list all csv file names in the folder and subfolders

#generate a 3 column dataframe with specific column names
Photo.R <- data.frame(matrix(NA, nrow=length(file.names)*2, ncol=5))
colnames(Photo.R) <- c("Fragment.ID","Intercept", "umol.L.sec","Temp.C","PR")

#Load Sample Info T1
Sample.Info <- read.csv(file="data/2017/Colony_Respirometry/Colony_Sample_Info_T1.csv", header=T) #read in volume and sample.info data
Sample.Info$Vol.L <- Sample.Info$Chamber.Vol.mL/1000 #calculate volume

#PHOTOSYNTHESIS
path.p<-"data/2017/Colony_Respirometry/T1_clay/" #the location of all your respirometry files
file.names<-list.files(path = path.p, pattern = "csv$") #list all csv file names in the folder
Photo.R <- data.frame(matrix(NA, nrow=length(file.names)*2, ncol=3)) #generate a 3 column dataframe with specific column names
colnames(Photo.R) <- c("Fragment.ID","Intercept", "umol.L.sec")
file.names
#Add names for photosynthesis or respiration for for loop
PR<-c('Photo','Resp')

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Photo.Data1 <-read.table(file.path(path.p,file.names[i]), skip = 1, header=T, sep=",", na.string="NA", fill = TRUE, as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Photo.Data1  <- Photo.Data1[,c(2,9,16)] #subset columns of interest
  Photo.Data1$Time <- as.POSIXct(Photo.Data1$Time,format="%H:%M:%S", tz = "") #convert time from character to time
  brk <- which(diff(Photo.Data1$Time) > 30) #look for breaks in time of 30 seconds or more
  Photo <- subset(Photo.Data1, as.numeric(rownames(Photo.Data1)) < brk)  #subset by break in time stamp keeping everything before break
  Resp <- subset(Photo.Data1, as.numeric(rownames(Photo.Data1)) > brk) #subset by break in time stamp keeping everything before break
  lt.levs <- list(Photo, Resp) #list levels of segmentation

  for(j in 1:length(lt.levs)){
    Photo.Data <- as.data.frame(lt.levs[j])
    n<-dim(Photo.Data )[1] #identify length of data
    Photo.Data <-Photo.Data [120:(n-3),] #start at data point ~2 minute in to avoid excess noise from start of run and remove last 3 lines containing text
    n<-dim(Photo.Data )[1] #list length of trimmed data
    Photo.Data $sec <- 1:n #set seconds by one from start to finish of run

    #Save plot prior to and after data thinning to make sure thinning is not too extreme
    rename <- sub("_.*", "", file.names[i])
    pdf(paste0("output/Photo_Resp_Output/T1",rename,"_",j,"thinning.pdf"))
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

    #save original unthinned data
    Photo.Data.orig <- Photo.Data
    Photo.Data   <-  thinData(Photo.Data , by=20)$newData1 #thin data by every 20 points for all the O2 values
    Photo.Data$sec <- as.numeric(rownames(Photo.Data )) #maintain numeric values for time

    plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot thinned data
    usr  <-  par('usr') #plotting graphics using 'usr'
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA) #giving specific coordinates for plot using 'usr'
    whiteGrid()
    box()
    points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
    axis(1)
    axis(2, las=1)
    dev.off()
    ##Olito et al. 2017: It is running a bootstrapping technique and calculating the rate based on density
    Regs  <-  rankLocReg(xall=Photo.Data $sec, yall=Photo.Data $Value, alpha=0.2,
                         method="pc", verbose=TRUE)
    pdf(paste0("output/Photo_Resp_Output/T1",rename,"_",j,"regression.pdf"))
    plot(Regs)
    dev.off()

    s <- seq(0,nrow(Photo.R),length(lt.levs)) #to order the file output sequence in correct order in data frame
    Photo.R[j+s[i],2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
    Photo.R[j+s[i],1] <- rename #stores the file name in the Date column
    #Photo.R[j+s[i],1] <- paste0(rename,"_",j) #stores the file name in the Date column
    Photo.R[j+s[i],4] <- mean(Photo.Data$Temp, na.rm=T)  #stores the Temperature in the Temp.C column
    Photo.R[j+s[i],5] <- PR[j] #stores whether it is photosynthesis or respiration
  }
}

write.csv(Photo.R, 'output/Photo_Resp_Output/T1_Photo.R.csv')

Photo.R <- read.csv('output/Photo_Resp_Output/T1_Photo.R.csv')

#Split up the photostynthesis and respiration data into two dataframes
PHO <- Photo.R[Photo.R$V5=='Photo', ]
RES <- Photo.R[Photo.R$V5=='Resp', ]

#Load Sample Info
Sample.Info <- read.csv(file="data/2017/Colony_Respirometry/Colony_Sample_Info_T1.csv", header=T) #read in volume and sample.info data
#Convert sample volume to mL
Sample.Info$Vol.L <- Sample.Info$Chamber.Vol.mL/1000 #calculate volume

#Merge data with sample info
Resp <- merge(RES,Sample.Info, by="Fragment.ID" )
Photo <- merge(PHO,Sample.Info, by="Fragment.ID")

#Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
Resp$umol.sec <- Resp$umol.L.sec*Sample.Info$Vol.L
Photo$umol.sec <- Photo$umol.L.sec*Sample.Info$Vol.L

## BLANKS HAVE TO BE SPECIFIC TO RESPONSE VARIABLE (I.E., PHOTO OR RESP) AND TEMP
photo.blnk <- aggregate(umol.sec ~ Origin, data=Photo, mean)
photo.Blanks <- subset(photo.blnk, Origin == "Blank")
Photo$Blank <- photo.Blanks[1,2]

#Calculate resp blank rate
resp.blnk <- aggregate(umol.sec ~ Origin, data=Resp, mean)
resp.Blanks <- subset(resp.blnk, Origin == "Blank")
Resp$Blank <- resp.Blanks[1,2]

#Account for blank rate Subtract Blank by the temperature blank
Resp$umol.sec.corr <- Resp$umol.sec-Resp$Blank
Photo$umol.sec.corr <- Photo$umol.sec-Photo$Blank

#normalize to surface area and h-1
Resp$umol.cm2.hr <- (Resp$umol.sec.corr*3600)/Resp$Surf.Area.cm2
Photo$umol.cm2.hr <- (Photo$umol.sec.corr*3600)/Photo$Surf.Area.cm2

##T1 Results
#remove blanks from dataset
Photo <- subset(Photo, Origin!= "Blank")
Resp <- subset(Resp, Origin!= "Blank")

write.csv(Photo, file="output/Photo_Resp_Output/Photosynthesis.rates.T1.csv")
write.csv(Resp, file="output/Photo_Resp_Output/Respiration.rates.T1.csv")

#calculate gross photosynthesis Pnet -- Rdark

Photo <- read.csv("output/Photo_Resp_Output/Photosynthesis.rates.T1.csv")
Resp <- read.csv("output/Photo_Resp_Output/Respiration.rates.T1.csv")

resp.data <- merge(Photo[,c(2,12,13,14,21)],Resp[,c(2,21)], by="Fragment.ID")

#rename the columns
names(resp.data)[names(resp.data) == "umol.cm2.hr.x"]<- "Pnet_umol.cm2.hr" 
names(resp.data)[names(resp.data) == "umol.cm2.hr.y"] <- "Rdark_umol.cm2.hr"

#Pnet plus resp (if positive) is pGross
resp.data$Pgross_umol.cm2.hr <- resp.data$Pnet_umol.cm2.hr-resp.data$Rdark_umol.cm2.hr

write.csv(resp.data, file="output/Photo_Resp_Output/photo_resp.T1.csv")


#T2
#PHOTOSYNTHESIS
path.p<-"data/2017/Colony_Respirometry/T2" #the location of all your respirometry files
file.names<-list.files(path = path.p, pattern = "csv$") #list all csv file names in the folder
Photo.R <- data.frame(matrix(NA, nrow=length(file.names)*2, ncol=3)) #generate a 3 column dataframe with specific column names
colnames(Photo.R) <- c("Fragment.ID","Intercept", "umol.L.sec")

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Photo.Data1 <-read.table(file.path(path.p,file.names[i]), skip = 1, header=T, sep=",", na.string="NA", fill = TRUE, as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Photo.Data1  <- Photo.Data1[,c(2,9,16)] #subset columns of interest
  Photo.Data1$Time <- as.POSIXct(Photo.Data1$Time,format="%H:%M:%S", tz = "") #convert time from character to time
  brk <- which(diff(Photo.Data1$Time) > 30) #look for breaks in time of 30 seconds or more
  Photo <- subset(Photo.Data1, as.numeric(rownames(Photo.Data1)) < brk)  #subset by break in time stamp keeping everything before break
  Resp <- subset(Photo.Data1, as.numeric(rownames(Photo.Data1)) > brk) #subset by break in time stamp keeping everything before break
  lt.levs <- list(Photo, Resp) #list levels of segmentation
  
  for(j in 1:length(lt.levs)){    
    Photo.Data <- as.data.frame(lt.levs[j])
    n<-dim(Photo.Data )[1] #identify length of data
    Photo.Data <-Photo.Data [120:(n-3),] #start at data point ~2 minute in to avoid excess noise from start of run and remove last 3 lines containing text
    n<-dim(Photo.Data )[1] #list length of trimmed data
    Photo.Data $sec <- 1:n #set seconds by one from start to finish of run
    
    #Save plot prior to and after data thinning to make sure thinning is not too extreme
    rename <- sub("_.*", "", file.names[i])
    pdf(paste0("output/Colony_Photo_Resp_Output/T2",rename,"_",j,"thinning.pdf"))
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
    
    #save original unthinned data
    Photo.Data.orig <- Photo.Data
    Photo.Data   <-  thinData(Photo.Data , by=20)$newData1 #thin data by every 20 points for all the O2 values
    Photo.Data$sec <- as.numeric(rownames(Photo.Data )) #maintain numeric values for time
    
    plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot thinned data
    usr  <-  par('usr') #plotting graphics using 'usr'
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA) #giving specific coordinates for plot using 'usr'
    whiteGrid()
    box()
    points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
    axis(1)
    axis(2, las=1)
    dev.off()
    ##Olito et al. 2017: It is running a bootstrapping technique and calculating the rate based on density
    Regs  <-  rankLocReg(xall=Photo.Data $sec, yall=Photo.Data $Value, alpha=0.2, 
                         method="pc", verbose=TRUE) 
    pdf(paste0("output/Colony_Photo_Resp_Output/T2",rename,"_",j,"regression.pdf"))
    plot(Regs)
    dev.off()
    
    s <- seq(0,nrow(Photo.R),length(lt.levs)) #to order the file output sequence in correct order in data frame
    Photo.R[j+s[i],2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
    Photo.R[j+s[i],1] <- rename #stores the file name in the Date column
    #Photo.R[j+s[i],1] <- paste0(rename,"_",j) #stores the file name in the Date column
    Photo.R[j+s[i],4] <- mean(Photo.Data$Temp, na.rm=T)  #stores the Temperature in the Temp.C column
    Photo.R[j+s[i],5] <- PR[j] #stores whether it is photosynthesis or respiration
  }
}

write.csv(Photo.R, 'output/Photo_Resp_Output/T2_Photo.R.csv')

Photo.R <- read.csv('output/Photo_Resp_Output/T2_Photo.R.csv')

#Split up the photostynthesis and respiration data into two dataframes
PHO <- Photo.R[Photo.R$V5=='Photo', ]
RES <- Photo.R[Photo.R$V5=='Resp', ]

#Load Sample Info T2
Sample.Info2 <- read.csv(file="Colony_Sample_Info_T2.csv", header=T) #read in volume and sample.info data
Sample.Info2$Vol.L <- Sample.Info2$Chamber.Vol.mL/1000 #calculate volume
Sample.Info2

#Merge with sample info
Resp2 <- merge(RES,Sample.Info2, by="Fragment.ID")
Photo2 <- merge(PHO,Sample.Info2, by="Fragment.ID")

#split sample info between photo and resp
Sample.Info.P2<-Sample.Info2[Sample.Info2$Light_Dark=='Light',]
Sample.Info.R2<-Sample.Info2[Sample.Info2$Light_Dark=='Dark',]

#Account for chamber volume umol s-1
Resp2$umol.sec <- Resp2$umol.L.sec*Resp2$Vol.L
Photo2$umol.sec <- Photo2$umol.L.sec*Photo2$Vol.L

##Getting Blank means and attaching to data
photo.blnk <- aggregate(umol.sec ~ Origin*Treatment, data=Photo2, mean)
photo.Blanks <- subset(photo.blnk, Origin == "Blank") #values for blank
Photo2$Blank <- photo.Blanks$umol.sec[match(Photo2$Treatment,photo.Blanks$Treatment)] #matching the right treatment blank to samples

#Calculate resp blank rate
resp.blnk <- aggregate(umol.sec ~ Origin*Treatment, data=Resp2, mean)
resp.Blanks <- subset(resp.blnk, Origin == "Blank")
Resp2$Blank <- resp.Blanks$umol.sec[match(Resp2$Treatment,resp.Blanks$Treatment)]

#Account for blank rate Subtract Blank by the temperature blank
Resp2$umol.sec.corr <- Resp2$umol.sec-Resp2$Blank
Photo2$umol.sec.corr <- Photo2$umol.sec-Photo2$Blank

#normalize to surface area and h-1
Resp2$umol.cm2.hr <- (Resp2$umol.sec.corr*3600)/Resp2$Surf.Area.cm2
Photo2$umol.cm2.hr <- (Photo2$umol.sec.corr*3600)/Photo2$Surf.Area.cm2

#T2 Results

Photo2 <- read.csv("output/Photo_Resp_Output/Photosynthesis.rates.T2.csv")
Resp2 <- read.csv("output/Photo_Resp_Output/Respiration.rates.T2.csv")
#remove blanks
Photo2 <- subset(Photo2, Origin!= "Blank")
Resp2 <- subset(Resp2, Origin!= "Blank")

#calculate gross photosynthesis Pnet -- Rdark

resp.data <- merge(Photo2[,c(2,11,12,13,20)],Resp2[,c(2,20)], by="Fragment.ID")

#rename the column
names(resp.data)[names(resp.data) == "umol.cm2.hr.x"]<- "Pnet_umol.cm2.hr" 
names(resp.data)[names(resp.data) == "umol.cm2.hr.y"] <- "Rdark_umol.cm2.hr"

#Pnet plus resp (if positive) is pGross
resp.data$Pgross_umol.cm2.hr <- resp.data$Pnet_umol.cm2.hr-resp.data$Rdark_umol.cm2.hr
write.csv(resp.data, file="output/Photo_Resp_Output/photo_resp.T2.csv")

#### PLOTTING ####
# Plotting T1 base graph
resp.2017.T1 <- read.csv("output/Photo_Resp_Output/photo_resp.T1.csv")
resp.2017.T2 <- read.csv("output/Photo_Resp_Output/photo_resp.T2.csv")

resp.2017.T1$timepoint <- "Pre-Treatment"
resp.2017.T2$timepoint <- "Post-Treatment"

Resp.2017.all <- rbind(resp.2017.T1, resp.2017.T2)
Resp.2017.all$abs.resp <- -(Resp.2017.all$Rdark_umol.cm2.hr)

#Gross Photosynthesis T1
mean.A2017.gross.photo.T1 <- summarySE(Resp.2017.all, measurevar="Pgross_umol.cm2.hr", groupvars=c("timepoint", "Origin"))
mean.A2017.gross.photo.T1$timepoint <- factor(mean.A2017.gross.photo.T1$timepoint, levels = c("Pre-Treatment","Post-Treatment"))

# Plotting
pd <- position_dodge(0.3) # moves object .05 to the left and right
#legend.title <- "Treatment"
PGross2017Adult.T1 <- ggplot(mean.A2017.gross.photo.T1, aes(x=timepoint, y=Pgross_umol.cm2.hr, group=Origin)) + 
  #  geom_line(position=pd, color="black", aes(linetype=reef.zone), size = 2)+
  geom_errorbar(aes(ymin=Pgross_umol.cm2.hr-se, ymax=Pgross_umol.cm2.hr+se), width=.1, size = 1, position=pd, color="black") + #Error bars
#  ylim(0,5.0)+
  xlab("Timepoint") + ylab(expression("Gross Photosynthetic Rate " (mu*mol ~ cm^{-2} ~ h^{-1}))) + #Axis titles
  geom_point(aes(shape=Origin), fill = "black", size=14, position=pd, color = "black")+
  scale_shape_manual(values=c(21,24),
                     name = "Reef Zone")+
  #  scale_fill_manual(values=c("#FFFFFF", "#999999"),
  #                    name = "Treatment")+ #colour modification
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  annotate("rect", xmin=0, xmax=1.5, ymin = 0, ymax=1.2, alpha = 0.2) + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/A2017.PGross.T1.pdf", PGross2017Adult.T1, width = 11, height = 11, units = c("in"))

# Pgross T1 Statistics

A2017.respall.T1 <- Resp.2017.all %>%
  filter(timepoint == "Pre-Treatment")

A2017.Pgross.T1.anova <- lm(Pgross_umol.cm2.hr~Origin, data = A2017.respall.T1)
qqnorm(resid(A2017.Pgross.T1.anova))
qqline(resid(A2017.Pgross.T1.anova))

boxplot(resid(A2017.Pgross.T1.anova)~A2017.respall.T1$Origin) # not equal


wilcox.test(Pgross_umol.cm2.hr~Origin, data = A2017.respall.T1)

capture.output(wilcox.test(Pgross_umol.cm2.hr~Origin, data = A2017.respall.T1), file = "output/Statistics/A2017.Pgross.T1.csv")



# Respiration T1 

mean.A2017.resp.T1 <- summarySE(Resp.2017.all, measurevar="abs.resp", groupvars=c("timepoint", "Origin"))
mean.A2017.resp.T1$timepoint <- factor(mean.A2017.resp.T1$timepoint, levels = c("Pre-Treatment","Post-Treatment"))

# Plotting
pd <- position_dodge(0.3) # moves object .05 to the left and right
#legend.title <- "Treatment"
Resp2017Adult.T1 <- ggplot(mean.A2017.resp.T1, aes(x=timepoint, y=abs.resp, group=Origin)) + 
  #  geom_line(position=pd, color="black", aes(linetype=reef.zone), size = 2)+
  geom_errorbar(aes(ymin=abs.resp-se, ymax=abs.resp+se), width=.1, size = 1, position=pd, color="black") + #Error bars
  #  ylim(0,5.0)+
  xlab("Timepoint") + ylab(expression("Respiration Rate " (mu*mol ~ cm^{-2} ~ h^{-1}))) + #Axis titles
  geom_point(aes(shape=Origin), fill = "black", size=14, position=pd, color = "black")+
  scale_shape_manual(values=c(21,24),
                     name = "Reef Zone")+
  #  scale_fill_manual(values=c("#FFFFFF", "#999999"),
  #                    name = "Treatment")+ #colour modification
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  annotate("rect", xmin=0, xmax=1.5, ymin = 0.2, ymax=0.7, alpha = 0.2) + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/A2017.Resp.T1.pdf", Resp2017Adult.T1, width = 11, height = 11, units = c("in"))

# T1 Resp Statsistics 

A2017.respall.T1 <- Resp.2017.all %>%
  filter(timepoint == "Pre-Treatment")

A2017.Resp.T1.anova <- lm(abs.resp~Origin, data = A2017.respall.T1)
qqnorm(resid(A2017.Resp.T1.anova))
qqline(resid(A2017.Resp.T1.anova))

boxplot(resid(A2017.Resp.T1.anova)~A2017.respall.T1$Origin) # not equal

wilcox.test(abs.resp~Origin, data = A2017.respall.T1)

capture.output(wilcox.test(abs.resp~Origin, data = A2017.respall.T1), file = "output/Statistics/A2017.Resp.T1.csv")


# Plotting T2 

A2017.resp.all.T2 <- Resp.2017.all %>%
  filter(timepoint == "Post-Treatment")


# Gross photosynthesis 
mean.A2017.Pgross.T2 <- summarySE(A2017.resp.all.T2, measurevar="Pgross_umol.cm2.hr", groupvars=c("Treatment", "Origin"))

# Plotting
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Treatment"
PGross2017Adult.T2 <- ggplot(mean.A2017.Pgross.T2, aes(x=Treatment, y=Pgross_umol.cm2.hr, group=Origin)) + 
  geom_line(position=pd, color="black", aes(linetype=Origin), size = 2)+
  geom_errorbar(aes(ymin=Pgross_umol.cm2.hr-se, ymax=Pgross_umol.cm2.hr+se), width=.1, size = 1, position=pd, color="black") + #Error bars
  ylim(0,1)+
  xlab("Timepoint") + ylab(expression("Gross Photosynthetic Rate " (mu*mol ~ cm^{-2} ~ h^{-1})))+ #Axis titles
  geom_point(aes(fill=Treatment, shape=Origin), size=14, position=pd, color = "black")+
  scale_shape_manual(values=c(21,24),
                     name = "Reef Zone")+
  scale_fill_manual(values=c("#FFFFFF", "#999999"),
                    name = "Treatment")+ #colour modification
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  annotate("rect", xmin=0, xmax=0.75, ymin = 0, ymax=1.2, alpha = 0.2) + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/A2017.PGross.T2.pdf", PGross2017Adult.T2, width = 11, height = 11, units = c("in"))


# Respiration Rate 
mean.A2017.Resp.T2 <- summarySE(A2017.resp.all.T2, measurevar="abs.resp", groupvars=c("Treatment", "Origin"))

# Plotting
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Treatment"
Resp2017Adult.T2 <- ggplot(mean.A2017.Resp.T2, aes(x=Treatment, y=abs.resp, group=Origin)) + 
  geom_line(position=pd, color="black", aes(linetype=Origin), size = 2)+
  geom_errorbar(aes(ymin=abs.resp-se, ymax=abs.resp+se), width=.1, size = 1, position=pd, color="black") + #Error bars
  ylim(0,1)+
  xlab("Timepoint") +  ylab(expression("Respiration Rate " (mu*mol ~ cm^{-2} ~ h^{-1})))+ #Axis titles
  geom_point(aes(fill=Treatment, shape=Origin), size=14, position=pd, color = "black")+
  scale_shape_manual(values=c(21,24),
                     name = "Reef Zone")+
  scale_fill_manual(values=c("#FFFFFF", "#999999"),
                    name = "Treatment")+ #colour modification
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  annotate("rect", xmin=0, xmax=0.75, ymin = 0.2, ymax=0.7, alpha = 0.2) + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(colour = "black", size=1), legend.position = "none") +
  theme(axis.text = element_text(size = 30, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

ggsave(file = "output/Graphs/A2017.Resp.T2.pdf", Resp2017Adult.T2, width = 11, height = 11, units = c("in"))


## Statsitics (time 2)
## Statistics

#Pnet
Pnet2017anova <- lm(Pnet_umol.cm2.hr~Origin*Treatment, data = A2017.resp.all.T2)
qqnorm(resid(Pnet2017anova))
qqline(resid(Pnet2017anova)) 

boxplot(resid(Pnet2017anova)~A2017.resp.all.T2$Origin)
boxplot(resid(Pnet2017anova)~A2017.resp.all.T2$Treatment) 

anova(Pnet2017anova)

capture.output(anova(Pnet2017anova), file = "output/Statistics/A2017.Pnet.csv")

#Pgross
Pgross2017anova <- lm(Pgross_umol.cm2.hr~Origin*Treatment, data = A2017.resp.all.T2)
qqnorm(resid(Pgross2017anova))
qqline(resid(Pgross2017anova)) 

boxplot(resid(Pgross2017anova)~A2017.resp.all.T2$Origin)
boxplot(resid(Pgross2017anova)~A2017.resp.all.T2$Treatment) 

anova(Pgross2017anova)

capture.output(anova(Pgross2017anova), file = "output/Statistics/A2017.Pgross.csv")

#Rdark
Rdark2017anova <- lm(Rdark_umol.cm2.hr~Origin*Treatment, data = A2017.resp.all.T2)
qqnorm(resid(Rdark2017anova))
qqline(resid(Rdark2017anova)) 

boxplot(resid(Rdark2017anova)~A2017.resp.all.T2$Origin)
boxplot(resid(Rdark2017anova)~A2017.resp.all.T2$Treatment) 

anova(Rdark2017anova)

capture.output(anova(Rdark2017anova), file = "output/Statistics/A2017.Rdark.csv")
