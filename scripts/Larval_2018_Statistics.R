#Title: Thermal Transplant 2017-2018 Statsitics (Larvae 2018)
#Author: KH Wong
#Date Last Modified: 20200711
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
library(lme4)
library(jtools)
library(car)

#### Complex experimental design with 3 fixed factors (Origin, treatment, transplant site) and potential random factors (Date and Coral ID)
#### Need to compare all 3 models, with different combinations of random factors

## Model comparison functionn 
#two *nested* models as input, in this case obj1 is the model without the predictor variable and obj2 is the model with the predictor variable 
lrt <- function (obj1, obj2) {
  L0 <- logLik(obj1)
  L1 <- logLik(obj2)
  L01 <- as.vector(- 2 * (L0 - L1))
  df <- attr(L1, "df") - attr(L0, "df")
  list(L01 = L01, df = df,
       "p-value" = pchisq(L01, df, lower.tail = FALSE))
}

#### CELL DENSITY #### 

L2018.zoox.mean.col.size <- read.csv("data/2018/Zoox/L2018_Zoox_calc.csv")


### orthogonal random effect design 

#Three way model with date and coral ID as random factor
Zoox2018Larvae.2018.lmer <- lmer(Cells.x3.mm3~Origin*Treatment.y*Transplant.Site + (1|Date)  +  (1|coral.id.x), data = L2018.zoox.mean.col.size, REML=FALSE)

qqnorm(resid(Zoox2018Larvae.2018.lmer)) # Normality
qqline(resid(Zoox2018Larvae.2018.lmer)) # Normal

boxplot(resid(Zoox2018Larvae.2018.lmer)~L2018.zoox.mean.col.size$Origin) # Variance 
boxplot(resid(Zoox2018Larvae.2018.lmer)~L2018.zoox.mean.col.size$Treatment.y)
boxplot(resid(Zoox2018Larvae.2018.lmer)~L2018.zoox.mean.col.size$Transplant.Site)

summary(Zoox2018Larvae.2018.lmer)

#Three way model with date as random factor
Zoox2018Larvae.2018.lmer2 <- lmer(Cells.x3.mm3~Origin*Treatment.y*Transplant.Site + (1|Date) , data = L2018.zoox.mean.col.size, REML=FALSE)
qqnorm(resid(Zoox2018Larvae.2018.lmer2)) # Normality
qqline(resid(Zoox2018Larvae.2018.lmer2)) # Normal

boxplot(resid(Zoox2018Larvae.2018.lmer2)~L2018.zoox.mean.col.size$Origin) # Variance 
boxplot(resid(Zoox2018Larvae.2018.lmer2)~L2018.zoox.mean.col.size$Treatment.y)
boxplot(resid(Zoox2018Larvae.2018.lmer2)~L2018.zoox.mean.col.size$Transplant.Site)

summary(Zoox2018Larvae.2018.lmer2)


#Three way model with coral.id as random factor
Zoox2018Larvae.2018.lmer3 <- lmer(Cells.x3.mm3~Origin*Treatment.y*Transplant.Site + (1|coral.id.x) , data = L2018.zoox.mean.col.size, REML=FALSE)
qqnorm(resid(Zoox2018Larvae.2018.lmer3)) # Normality
qqline(resid(Zoox2018Larvae.2018.lmer3)) # Normal

boxplot(resid(Zoox2018Larvae.2018.lmer3)~L2018.zoox.mean.col.size$Origin) # Variance 
boxplot(resid(Zoox2018Larvae.2018.lmer3)~L2018.zoox.mean.col.size$Treatment.y)
boxplot(resid(Zoox2018Larvae.2018.lmer3)~L2018.zoox.mean.col.size$Transplant.Site)

summary(Zoox2018Larvae.2018.lmer3)

## Model Comparisons
lrt(Zoox2018Larvae.2018.lmer2, Zoox2018Larvae.2018.lmer) #Full model is best
lrt(Zoox2018Larvae.2018.lmer3, Zoox2018Larvae.2018.lmer) #Full model is best

#MODEL SELECTION = FULL 
capture.output(Anova(Zoox2018Larvae.2018.lmer), file = "output/Statistics/L.2018.zoox.lmer.csv")


#### VOLUME #### 

L2018.vol.mean.ColDay <- read.csv("data/2018/Larval.Size/L.Size.2018.COLDAY.csv")

### orthogonal random effect design 

#Three way model with date and coral ID as random factor
Vol2018Larvae.2018.lmer <- lmer(Volume~Origin*Treatment*Transplant.Site + (1|Date.x) +  (1|coral.id), data = L2018.vol.mean.ColDay, REML=FALSE)
qqnorm(resid(Vol2018Larvae.2018.lmer)) # Normality
qqline(resid(Vol2018Larvae.2018.lmer)) # Normal

boxplot(resid(Vol2018Larvae.2018.lmer)~L2018.vol.mean.ColDay$Origin) # Variance 
boxplot(resid(Vol2018Larvae.2018.lmer)~L2018.vol.mean.ColDay$Treatment)
boxplot(resid(Vol2018Larvae.2018.lmer)~L2018.vol.mean.ColDay$Transplant.Site)

summary(Vol2018Larvae.2018.lmer)

#Three way model with date as random factor
Vol2018Larvae.2018.lmer2 <- lmer(Volume~Origin*Treatment*Transplant.Site + (1|Date.x), data = L2018.vol.mean.ColDay, REML=FALSE)
qqnorm(resid(Vol2018Larvae.2018.lmer2)) # Normality
qqline(resid(Vol2018Larvae.2018.lmer2)) # Normal

boxplot(resid(Vol2018Larvae.2018.lmer2)~L2018.vol.mean.ColDay$Origin) # Variance 
boxplot(resid(Vol2018Larvae.2018.lmer2)~L2018.vol.mean.ColDay$Treatment)
boxplot(resid(Vol2018Larvae.2018.lmer2)~L2018.vol.mean.ColDay$Transplant.Site)

summary(Vol2018Larvae.2018.lmer2)

#Three way model with coral.id as random factor
Vol2018Larvae.2018.lmer3 <- lmer(Volume~Origin*Treatment*Transplant.Site + (1|coral.id), data = L2018.vol.mean.ColDay, REML=FALSE)
qqnorm(resid(Vol2018Larvae.2018.lmer3)) # Normality
qqline(resid(Vol2018Larvae.2018.lmer3)) # Normal

boxplot(resid(Vol2018Larvae.2018.lmer3)~L2018.vol.mean.ColDay$Origin) # Variance 
boxplot(resid(Vol2018Larvae.2018.lmer3)~L2018.vol.mean.ColDay$Treatment)
boxplot(resid(Vol2018Larvae.2018.lmer3)~L2018.vol.mean.ColDay$Transplant.Site)

summary(Vol2018Larvae.2018.lmer3)

## Model Comparisons
lrt(Vol2018Larvae.2018.lmer2, Vol2018Larvae.2018.lmer) #Model II is best
lrt(Vol2018Larvae.2018.lmer3, Vol2018Larvae.2018.lmer) #Full model is best

#MODEL SELECTION = II 
capture.output(Anova(Vol2018Larvae.2018.lmer2), file = "output/Statistics/L.2018.vol.lmer.csv")


#### TOTAL PROTEIN #### 

L2018.TP.mean.size.col.day <- read.csv("data/2018/Protein/L2018.TP.Vol.Means.csv")

### orthogonal random effect design 

#Three way model with date and coral ID as random factor
TP2018Larvae.2018.lmer <- lmer(Conc.calcS.mg.mm3~Origin*Treatment.y*Transplant.Site + (1|Date)+ (1|coral.id.x), data = L2018.TP.mean.size.col.day, REML=FALSE)
qqnorm(resid(TP2018Larvae.2018.lmer)) # Normality
qqline(resid(TP2018Larvae.2018.lmer)) # Normal

boxplot(resid(TP2018Larvae.2018.lmer)~L2018.TP.mean.size.col.day$Origin) # Variance 
boxplot(resid(TP2018Larvae.2018.lmer)~L2018.TP.mean.size.col.day$Treatment)
boxplot(resid(TP2018Larvae.2018.lmer)~L2018.TP.mean.size.col.day$Transplant.Site)

summary(TP2018Larvae.2018.lmer)

#Three way model with date as random factor
TP2018Larvae.2018.lmer2 <- lmer(Conc.calcS.mg.mm3~Origin*Treatment.y*Transplant.Site + (1|Date), data = L2018.TP.mean.size.col.day, REML=FALSE)
qqnorm(resid(TP2018Larvae.2018.lmer2)) # Normality
qqline(resid(TP2018Larvae.2018.lmer2)) # Normal

boxplot(resid(TP2018Larvae.2018.lmer2)~L2018.TP.mean.size.col.day$Origin) # Variance 
boxplot(resid(TP2018Larvae.2018.lmer2)~L2018.TP.mean.size.col.day$Treatment)
boxplot(resid(TP2018Larvae.2018.lmer2)~L2018.TP.mean.size.col.day$Transplant.Site)

summary(TP2018Larvae.2018.lmer2)

#Three way model with coral.id as random factor
TP2018Larvae.2018.lmer3 <- lmer(Conc.calcS.mg.mm3~Origin*Treatment.y*Transplant.Site + (1|coral.id.x), data = L2018.TP.mean.size.col.day, REML=FALSE)
qqnorm(resid(TP2018Larvae.2018.lmer3)) # Normality
qqline(resid(TP2018Larvae.2018.lmer3)) # Normal

boxplot(resid(TP2018Larvae.2018.lmer3)~L2018.TP.mean.size.col.day$Origin) # Variance 
boxplot(resid(TP2018Larvae.2018.lmer3)~L2018.TP.mean.size.col.day$Treatment)
boxplot(resid(TP2018Larvae.2018.lmer3)~L2018.TP.mean.size.col.day$Transplant.Site)

summary(TP2018Larvae.2018.lmer3)

## Model Comparisons
lrt(TP2018Larvae.2018.lmer2, TP2018Larvae.2018.lmer) #Model II is best
lrt(TP2018Larvae.2018.lmer3, TP2018Larvae.2018.lmer) #Full model is best

#MODEL SELECTION = II 
capture.output(Anova(TP2018Larvae.2018.lmer2), file = "output/Statistics/L.2018.TP.lmer.csv")

#### CHLA/VOL #### 

L2018.Chla.col.day <- read.csv("data/2018/Chlorophyll/L2018.Chlavol.mean.COLDAY.csv")

### orthogonal random effect design 

#Three way model with date and coral ID as random factor
Chla2018Larvae.2018.lmer <- lmer(Chla.ng.mm3~Origin*Treatment*Transplant.Site + (1|Date) + (1|coral.id.x), data = L2018.Chla.col.day, REML=FALSE)
qqnorm(resid(Chla2018Larvae.2018.lmer)) # Normality
qqline(resid(Chla2018Larvae.2018.lmer)) # Normal

boxplot(resid(Chla2018Larvae.2018.lmer)~L2018.Chla.col.day$Origin) # Variance 
boxplot(resid(Chla2018Larvae.2018.lmer)~L2018.Chla.col.day$Treatment)
boxplot(resid(Chla2018Larvae.2018.lmer)~L2018.Chla.col.day$Transplant.Site)

summary(Chla2018Larvae.2018.lmer)

#Three way model with date as random factor
Chla2018Larvae.2018.lmer2 <- lmer(Chla.ng.mm3~Origin*Treatment*Transplant.Site + (1|Date), data = L2018.Chla.col.day, REML=FALSE)
qqnorm(resid(Chla2018Larvae.2018.lmer2)) # Normality
qqline(resid(Chla2018Larvae.2018.lmer2)) # Normal

boxplot(resid(Chla2018Larvae.2018.lmer2)~L2018.Chla.col.day$Origin) # Variance 
boxplot(resid(Chla2018Larvae.2018.lmer2)~L2018.Chla.col.day$Treatment)
boxplot(resid(Chla2018Larvae.2018.lmer2)~L2018.Chla.col.day$Transplant.Site)

summary(Chla2018Larvae.2018.lmer2)

#Three way model with coral.id as random factor
Chla2018Larvae.2018.lmer3 <- lmer(Chla.ng.mm3~Origin*Treatment*Transplant.Site + (1|coral.id.x), data = L2018.Chla.col.day, REML=FALSE)
qqnorm(resid(Chla2018Larvae.2018.lmer3)) # Normality
qqline(resid(Chla2018Larvae.2018.lmer3)) # Normal

boxplot(resid(Chla2018Larvae.2018.lmer3)~L2018.Chla.col.day$Origin) # Variance 
boxplot(resid(Chla2018Larvae.2018.lmer3)~L2018.Chla.col.day$Treatment)
boxplot(resid(Chla2018Larvae.2018.lmer3)~L2018.Chla.col.day$Transplant.Site)

summary(Chla2018Larvae.2018.lmer3)

## Model Comparisons
lrt(Chla2018Larvae.2018.lmer2, Chla2018Larvae.2018.lmer) #Full model is best
lrt(Chla2018Larvae.2018.lmer3, Chla2018Larvae.2018.lmer) #Full model is best

#MODEL SELECTION = III 
capture.output(Anova(Chla2018Larvae.2018.lmer3), file = "output/Statistics/L.2018.Chlavol.lmer.csv")


#### CHLA/VOL #### 

L2018.Chlacell.col.day <- read.csv("data/2018/Chlorophyll/L2018.Chlacell.col.day.csv")

### orthogonal random effect design 

#Three way model with date and coral ID as random factor
Chlacell2018Larvae.2018.lmer <- lmer(chla.cell~Origin*Treatment*Transplant.Site + (1|Date) + (1|coral.id.x), data = L2018.Chlacell.col.day, REML=FALSE)
qqnorm(resid(Chlacell2018Larvae.2018.lmer)) # Normality
qqline(resid(Chlacell2018Larvae.2018.lmer)) # Normal

boxplot(resid(Chlacell2018Larvae.2018.lmer)~L2018.Chlacell.col.day$Origin) # Variance 
boxplot(resid(Chlacell2018Larvae.2018.lmer)~L2018.Chlacell.col.day$Treatment)
boxplot(resid(Chlacell2018Larvae.2018.lmer)~L2018.Chlacell.col.day$Transplant.Site)

summary(Chlacell2018Larvae.2018.lmer)

#Three way model with date as random factor
Chlacell2018Larvae.2018.lmer2 <- lmer(chla.cell~Origin*Treatment*Transplant.Site + (1|Date), data = L2018.Chlacell.col.day, REML=FALSE)
qqnorm(resid(Chlacell2018Larvae.2018.lmer2)) # Normality
qqline(resid(Chlacell2018Larvae.2018.lmer2)) # Normal

boxplot(resid(Chlacell2018Larvae.2018.lmer2)~L2018.Chlacell.col.day$Origin) # Variance 
boxplot(resid(Chlacell2018Larvae.2018.lmer2)~L2018.Chlacell.col.day$Treatment)
boxplot(resid(Chlacell2018Larvae.2018.lmer2)~L2018.Chlacell.col.day$Transplant.Site)

summary(Chlacell2018Larvae.2018.lmer2)

#Three way model with coral.id as random factor
Chlacell2018Larvae.2018.lmer3 <- lmer(chla.cell~Origin*Treatment*Transplant.Site + (1|coral.id.x), data = L2018.Chlacell.col.day, REML=FALSE)
qqnorm(resid(Chlacell2018Larvae.2018.lmer3)) # Normality
qqline(resid(Chlacell2018Larvae.2018.lmer3)) # Normal

boxplot(resid(Chlacell2018Larvae.2018.lmer3)~L2018.Chlacell.col.day$Origin) # Variance 
boxplot(resid(Chlacell2018Larvae.2018.lmer3)~L2018.Chlacell.col.day$Treatment)
boxplot(resid(Chlacell2018Larvae.2018.lmer3)~L2018.Chlacell.col.day$Transplant.Site)

summary(Chlacell2018Larvae.2018.lmer3)


## Model Comparisons
lrt(Chlacell2018Larvae.2018.lmer2, Chlacell2018Larvae.2018.lmer) #Model II is best
lrt(Chlacell2018Larvae.2018.lmer3, Chlacell2018Larvae.2018.lmer) #Model III is best
lrt(Chlacell2018Larvae.2018.lmer3, Chlacell2018Larvae.2018.lmer2) #Model II is best

#MODEL SELECTION = III 
capture.output(Anova(Chlacell2018Larvae.2018.lmer2), file = "output/Statistics/L.2018.Chlacell.lmer.csv")
