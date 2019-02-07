#Title: Photosynthesis and Respiration Calculations for 2017-2018 Transplant Experiment 
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

#PI.curve.data<-read.csv('~/Dropbox/BIOS_Summer_2018/20172018_Transplant/RAnalysis/Data/Transp.PI.Curve.csv', header=T, sep = ",")
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




##### PHOTOSYNTHESIS POST-TRANSPLANT (2018) #####

path.p<-"~/Dropbox/BIOS_Summer_2018/20172018_Transplant/RAnalysis/Data/OR_RESP/" #the location of all your respirometry files 

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
    pdf(paste0("~/Dropbox/BIOS_Summer_2018/20172018_Transplant/RAnalysis/Output/Photo_Resp_Output/",rename,"_",j,"thinning.pdf"))
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
    pdf(paste0("~/Dropbox/BIOS_Summer_2018/20172018_Transplant/RAnalysis/Output/Photo_Resp_Output/",rename,"_",j,"regression.pdf"))
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

write.csv(Photo.R, 'Photo.R.csv')

Photo.R <- read.csv(file="~/Dropbox/BIOS_Summer_2018/20172018_Transplant/RAnalysis/Data/Photo.R.csv", header=T) #read in volume and sample.info data

#Split up the photostynthesis and respiration data into two dataframes
PHO <- Photo.R[Photo.R$V5=='Photo', ]
RES <- Photo.R[Photo.R$V5=='Resp', ]

#Removing _1 or _2
PHO$Fragment.ID <- PHO$Fragment.ID %>% str_replace("_1","")
RES$Fragment.ID <- RES$Fragment.ID %>% str_replace("_2","")

#Load Sample Info
Sample.Info <- read.csv(file="~/Dropbox/BIOS_Summer_2018/20172018_Transplant/RAnalysis/Data/Sample_Info_Transp.csv", header=T) #read in volume and sample.info data
Sample.Info$Vol.L <- Sample.Info$Height*15.625*8.000*0.0163871 #calculate volume from height of water and dimensions of tank, multiplied by the conversion of cubic inch to liters
colnames(Sample.Info)[colnames(Sample.Info)=="coral.id"] <- "Fragment.ID" #changing column name to match resp data

coral.sa.calc <- read.csv('~/Dropbox/BIOS_Summer_2018/20172018_Transplant/RAnalysis/Data/coral.trans.calcsa.csv', header=T, stringsAsFactors=FALSE)
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
photo.blnk <- aggregate(umol.sec ~ Group, data=Photo, mean)
photo.Blanks <- subset(photo.blnk, Group == "Blank")
Photo$Blank <- photo.Blanks[1,2]

#Calculate resp blank rate
resp.blnk <- aggregate(umol.sec ~ Group, data=Resp, mean)
resp.Blanks <- subset(resp.blnk, Group == "Blank")
Resp$Blank <- resp.Blanks[1,2]

#Account for blank rate Subtract Blank by the temperature blank
Resp$umol.sec.corr <- Resp$umol.sec-Resp$Blank
Photo$umol.sec.corr <- Photo$umol.sec-Photo$Blank

#remove blanks from dataset
Photo <- subset(Photo, Group!= "Blank")
Resp <- subset(Resp, Group!= "Blank")

#normalize to surface area and h-1
Photo$umol.cm2.hr<-(Photo$umol.sec.corr*3600) / Photo$Surface_Area
Resp$umol.cm2.hr<-(Resp$umol.sec.corr*3600) / Resp$Surface_Area

## Results

write.csv(Photo, file="~/Dropbox/BIOS_Summer_2018/20172018_Transplant/RAnalysis/Output/Photo_Resp_Output/Photosynthesis.rates.csv")
write.csv(Resp, file="~/Dropbox/BIOS_Summer_2018/20172018_Transplant/RAnalysis/Output/Photo_Resp_Output/Respiration.rates.csv")

#calculate gross photosynthesis Pnet -- Rdark
resp.data <- merge(Photo[,c(1,11,12,13,14,26)],Resp[,c(1,26)], by="Fragment.ID")

#rename the columns
#names(resp.data)[names(resp.data) == "umol.cm2.hr.x"]<- "Pnet_umol.cm2.hr" 
#names(resp.data)[names(resp.data) == "umol.cm2.hr.y"] <- "Rdark_umol.cm2.hr"
colnames(resp.data)[6] <- "Pnet_umol.cm2.hr" #not necessary
colnames(resp.data) [7] <- "Rdark_umol.cm2.hr"

#Pnet plus resp (if positive) is pGross
resp.data$Pgross_umol.cm2.hr <- resp.data$Pnet_umol.cm2.hr - resp.data$Rdark_umol.cm2.hr

#adding treatment column
#resp.data <- merge(resp.data[,],Photo[,c(1,11)], by="Fragment.ID")

resp.data$Origin.Treatment <- paste(resp.data$Origin, resp.data$Treatment)

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
