#Title: Photosynthesis and Respiration Calculations for 2017-2018 Transplant Experiment 
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
library(stringr)

##### PHOTOSYNTHESIS POST-TRANSPLANT (2018) #####

path.p<-"data/2018/Colony_Respirometry/OR_RESP/" #the location of all your respirometry files 

# bring in the respiration files
file.names<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders
#basename above removes the subdirectory name from the file

#add file names that include the subdirectory name
file.names.full<-list.files(path = path.p, pattern = "csv$", recursive = TRUE) #list all csv file names in the folder and subfolders

#generate a 3 column dataframe with specific column names
Photo.R <- data.frame(matrix(NA, nrow=length(file.names)*2, ncol=5))
colnames(Photo.R) <- c("Fragment.ID","Intercept", "umol.L.sec","Temp.C","PR")

#PHOTOSYNTHESIS
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
    pdf(paste0("output/Photo_Resp_Output/",rename,"_",j,"thinning.pdf"))
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
    pdf(paste0("output/Photo_Resp_Output/",rename,"_",j,"regression.pdf"))
    plot(Regs)
    dev.off()
    
    s <- seq(0,nrow(Photo.R),length(lt.levs)) #to order the file output sequence in correct order in data frame
    Photo.R[j+s[i],2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
    Photo.R[j+s[i],1] <- rename #stores the file name in the Date column
    Photo.R[j+s[i],1] <- paste0(rename,"_",j) #stores the file name in the Date column
    Photo.R[j+s[i],4] <- mean(Photo.Data$Temp, na.rm=T)  #stores the Temperature in the Temp.C column
    Photo.R[j+s[i],5] <- PR[j] #stores whether it is photosynthesis or respiration
  }
}

write.csv(Photo.R, 'data/2018/Colony_Respirometry/Photo.R.csv')

Photo.R <- read.csv(file="data/2018/Colony_Respirometry/Photo.R.csv", header=T) #read in volume and sample.info data

#Split up the photostynthesis and respiration data into two dataframes
PHO <- Photo.R[Photo.R$V5=='Photo', ]
RES <- Photo.R[Photo.R$V5=='Resp', ]

#Removing _1 or _2
PHO$Fragment.ID <- PHO$Fragment.ID %>% str_replace("_1","")
RES$Fragment.ID <- RES$Fragment.ID %>% str_replace("_2","")

#Load Sample Info
Sample.Info <- read.csv(file="data/2018/Colony_Respirometry/Sample_Info_Transp.csv", header=T) #read in volume and sample.info data
Sample.Info$Vol.L <- Sample.Info$Height*15.625*8.000*0.0163871 #calculate volume from height of water and dimensions of tank, multiplied by the conversion of cubic inch to liters
colnames(Sample.Info)[colnames(Sample.Info)=="coral.id"] <- "Fragment.ID" #changing column name to match resp data

coral.sa.calc <- read.csv('data/2018/Surface.Area/colony.2018.calcsa.resp.csv', header=T, stringsAsFactors=FALSE)
colnames(coral.sa.calc)[colnames(coral.sa.calc)=="coral.id"] <- "Fragment.ID" #changing column name to match resp data
coral.sa.calc$Surface_Area <- as.numeric(coral.sa.calc$Surface_Area)
Sample.Info <- merge(Sample.Info, coral.sa.calc, by ="Fragment.ID")

#Merge data with sample info
Resp <- merge(RES,Sample.Info, by="Fragment.ID" )
Photo <- merge(PHO,Sample.Info, by="Fragment.ID")

#Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
Resp$umol.sec <- Resp$umol.L.sec*Sample.Info$Vol.L
Photo$umol.sec <- Photo$umol.L.sec*Sample.Info$Vol.L

## BLANKS HAVE TO BE SPECIFIC TO RESPONSE VARIABLE (I.E., PHOTO OR RESP) AND TEMP (ALL AT ONE TEMP)
photo.blnk <- aggregate(umol.sec ~ Type, data=Photo, mean)
photo.Blanks <- subset(photo.blnk, Type == "Blank")
Photo$Blank <- photo.Blanks[1,2]

#Calculate resp blank rate
resp.blnk <- aggregate(umol.sec ~ Type, data=Resp, mean)
resp.Blanks <- subset(resp.blnk, Type == "Blank")
Resp$Blank <- resp.Blanks[1,2]

#Account for blank rate Subtract Blank by the temperature blank
Resp$umol.sec.corr <- Resp$umol.sec-Resp$Blank
Photo$umol.sec.corr <- Photo$umol.sec-Photo$Blank

#remove blanks from dataset
Photo <- subset(Photo, Type!= "Blank")
Resp <- subset(Resp, Type!= "Blank")

#normalize to surface area and h-1
Photo$umol.cm2.hr<-(Photo$umol.sec.corr*3600) / Photo$Surface_Area
Resp$umol.cm2.hr<-(Resp$umol.sec.corr*3600) / Resp$Surface_Area

## Results

write.csv(Photo, file="output/Photo_Resp_Output/Photosynthesis.rates.2019.csv")
write.csv(Resp, file="output/Photo_Resp_Output/Respiration.rates2019.csv")

Photo <- read.csv("output/Photo_Resp_Output/Photosynthesis.rates.2019.csv")
Resp <- read.csv("output/Photo_Resp_Output/Respiration.rates2019.csv")

#calculate gross photosynthesis Pnet -- Rdark
resp.data <- merge(Photo[,c(1,10,11,12,13,25)],Resp[,c(1,25)], by="Fragment.ID")

#rename the columns
colnames(resp.data)[6] <- "Pnet_umol.cm2.hr" #not necessary
colnames(resp.data) [7] <- "Rdark_umol.cm2.hr"

#Pnet plus resp (if positive) is pGross
resp.data$Pgross_umol.cm2.hr <- resp.data$Pnet_umol.cm2.hr - resp.data$Rdark_umol.cm2.hr

resp.data$Origin.Treatment <- paste(resp.data$Origin, resp.data$Treatment)

write.csv(resp.data, file ="data/2018/Colony_Respirometry/All.resp.data2018.csv")

#Graphing

resp.data <- read.csv("data/2018/Colony_Respirometry/All.resp.data2018.csv")

### Gross Photosynthesis 

mean.PGross.A2018 <- summarySE(resp.data, measurevar="Pgross_umol.cm2.hr", groupvars=c("Origin", "Treatment", "Transplant.Site"))
mean.PGross.A2018$reef.treatment <- paste(mean.PGross.A2018$Origin, mean.PGross.A2018$Treatment)

#Plotting 
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Treatment"
PGross2018Adult <- ggplot(mean.PGross.A2018, aes(x=Transplant.Site, y=Pgross_umol.cm2.hr, group=reef.treatment)) + 
  geom_line(position=pd, aes(linetype=Origin, color = Treatment), size = 2)+
  geom_errorbar(aes(ymin=Pgross_umol.cm2.hr-se, ymax=Pgross_umol.cm2.hr+se), width=.1, position=pd, color="black") + #Error bars
  ylim(0.5,2.5)+
  xlab("Transplant Site") + ylab(expression("Gross Photosynthetic Rate " (mu*mol ~ cm^{-2} ~ h^{-1}))) + #Axis titles
  geom_point(aes(fill=Treatment, shape=Origin), size=14, position=pd, color = "black")+
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

PGross2018Adult
ggsave(file = "output/Graphs/A2018.PGross.pdf", PGross2018Adult,  width = 11, height = 11, units = c("in"))



### Respiration

resp.data$abs.resp <- -(resp.data$Rdark_umol.cm2.hr)
mean.Resp.A2018 <- summarySE(resp.data, measurevar="abs.resp", groupvars=c("Origin", "Treatment", "Transplant.Site"))
mean.Resp.A2018$reef.treatment <- paste(mean.Resp.A2018$Origin, mean.Resp.A2018$Treatment)

#Plotting 
pd <- position_dodge(0.1) # moves object .05 to the left and right
#legend.title <- "Treatment"
Resp2018Adult <- ggplot(mean.Resp.A2018, aes(x=Transplant.Site, y=abs.resp, group=reef.treatment)) + 
  geom_line(position=pd, aes(linetype=Origin, color = Treatment), size = 2)+
  geom_errorbar(aes(ymin=abs.resp-se, ymax=abs.resp+se), width=.1, position=pd, color="black") + #Error bars
  ylim(0.4,1.5)+
  xlab("Transplant Site") + ylab(expression("Respiration Rate " (mu*mol ~ cm^{-2} ~ h^{-1}))) + #Axis titles
  geom_point(aes(fill=Treatment, shape=Origin), size=14, position=pd, color = "black")+
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

Resp2018Adult
ggsave(file = "output/Graphs/A2018.Resp.pdf", Resp2018Adult,  width = 11, height = 11, units = c("in"))


### STATISTICS ###
#Pgross
Pgross2018anova <- lm(Pgross_umol.cm2.hr~Origin*Treatment*Transplant.Site, data = resp.data)
qqnorm(resid(Pgross2018anova))
qqline(resid(Pgross2018anova)) 

boxplot(resid(Pgross2018anova)~resp.data$Origin)
boxplot(resid(Pgross2018anova)~resp.data$Treatment) 
boxplot(resid(Pgross2018anova)~resp.data$Transplant.Site) 

anova(Pgross2018anova)

capture.output(anova(Pgross2018anova), file = "output/Statistics/A2018.Pgross.csv")

#Rdark
Rdark2018anova <- lm(Rdark_umol.cm2.hr~Origin*Treatment*Transplant.Site, data = resp.data)
qqnorm(resid(Rdark2018anova))
qqline(resid(Rdark2018anova)) 

boxplot(resid(Rdark2018anova)~resp.data$Origin)
boxplot(resid(Rdark2018anova)~resp.data$Treatment) 
boxplot(resid(Rdark2018anova)~resp.data$Transplant.Site) 

anova(Rdark2018anova)

capture.output(anova(Rdark2018anova), file = "output/Statistics/A2018.Rdark.csv")
