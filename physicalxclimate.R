#physical traits with proxy and fixed climatic variables 
rm(list=ls()) # clears work space

###install packages-------------------------------------------------------------
library(ggplot2)
library(MASS)
library(reshape2)
library(tidyverse)
library(dplyr)
library(randomForest)
library(glmmTMB) #mixed models
library(car) #mixed models summary tables
library(vegan) #multivariate stats
library(Boruta) #random forest models

#Read in data-------------------------------------------------------------------
library(readr)
c <- read_csv("climate_physical.csv")
View(c)

#How do management systems interact with broad abiotic conditions to shape fruit quality
#mgmt and phys traits-----------------------------------------------------------

m1<- glmmTMB(SSC ~ orchard.type + (1|site.code/onum/Tree), data=c)
summary(m1)
Anova(m1) #orchard.type p=0.294

m2<- glmmTMB(Firmness ~ orchard.type + (1|site.code/onum/Tree), data=c)
summary(m2)
Anova(m2)#orchard.type p=0.9851

m3<- glmmTMB(avgwgt ~ orchard.type + (1|site.code/onum/Tree), data=c)
summary(m3)
Anova(m3)#orchard.type p=0.8997

m4<- glmmTMB(maturity.index ~ orchard.type + (1|site.code/onum), data=c)
summary(m4)
Anova(m4)#orchard.type p=0.851

###within all of the orchards mgmt system alone doesn't have an effect on fruit physical traits 

###aggregating data for averages 
aggregate(SSC~orchard.type, data=c, FUN=mean)
# Conventional 10.88182
# Organic 11.78462

aggregate(Firmness~orchard.type, data=c, FUN=mean)
#Conventional 22.26473
#Organic 22.44308

aggregate(avgwgt~orchard.type, data=c, FUN=mean)
#1 Conventional 100.93273
#    Organic 98.68308

aggregate(maturity.index~orchard.type, data=c, FUN=mean)
# Conventional       3.781818
#    Organic       4.138462

###proxy climate and mgmt-------------------------------------------------------
###SSC###
s1 <- glmmTMB(SSC ~ Latitude*orchard.type + (1|site.code/onum/Tree), data=c)
summary(s1)
Anova(s1)
#Latitude p=1.906e-05 ***
#orchard.type p=0.4788    
#Latitude:orchard.type p=0.6393  

s2<- glmmTMB(SSC ~ elevation*orchard.type + (1|site.code/onum/Tree), data=c)
summary(s2)
Anova(s2)
#elevation p=0.1148
#orchard.type p=0.2849
#elevation:orchard.type p=0.6719

s3<- glmmTMB(SSC ~ Longitude*orchard.type + (1|site.code/onum/Tree), data=c)
summary(s3)
Anova(s3)
#Longitutde              0.5870  1     0.4436
#orchard.type            0.9519  1     0.3292
#Longitutde:orchard.type 0.0287  1     0.8656


###Firmness###
f1 <- glmmTMB(Firmness ~ Latitude*orchard.type + (1|site.code/onum/Tree), data=c)
summary(f1)
Anova(f1)
#Latitude p=4.112e-09 ***
#orchard.type p=0.7735    
#Latitude:orchard.type p=0.3366 

f2<- glmmTMB(Firmness ~ elevation*orchard.type + (1|site.code/onum/Tree), data=c)
summary(f2)
Anova(f2)
#elevation p=0.1808
#orchard.type p=0.9686
#elevation:orchard.type p=0.6368

f3<- glmmTMB(Firmness ~ Longitude*orchard.type + (1|site.code/onum/Tree), data=c)
summary(f3)
Anova(f3)

#Longitutde              0.2882  1     0.5914
#orchard.type            0.0030  1     0.9563
#Longitutde:orchard.type 1.8458  1     0.1743


###Average Weight (weight of 3 apples minus the bag divided by 3)###
w1 <- glmmTMB(avgwgt ~ Latitude*orchard.type + (1|site.code/onum/Tree), data=c)
summary(w1)
Anova(w1)
#Latitude p=0.570994   
#orchard.type p=0.716308   
#Latitude:orchard.type p=0.001643 ** 

w2<- glmmTMB(avgwgt ~ elevation*orchard.type + (1|site.code/onum/Tree), data=c)
summary(w2)
Anova(w2)
#elevation p=0.005707 **
#orchard.type p=0.465530   
#elevation:orchard.type p=0.087329 .

w3<- glmmTMB(avgwgt ~ Longitude*orchard.type + (1|site.code/onum/Tree), data=c)
summary(w3)
Anova(w3)
#Longitude p=0.3660
#orchard.type p=0.7945
#Longitude:orchard.type p=0.1693


#Maturity  
mt1 <- glmmTMB(maturity.index ~ Latitude*orchard.type + (1|site.code/onum/Tree), data=c)
summary(mt1)
Anova(mt1)
#Latitude p=0.0489 *
#orchard.type p=0.9454  
#Latitude:orchard.type p=0.3272 

mt2<- glmmTMB(maturity.index ~ elevation*orchard.type + (1|site.code/onum/Tree), data=c)
summary(mt2)
Anova(mt2)
#elevation p=0.5632
#orchard.type p=0.8615
#elevation:orchard.type p=0.4549

mt3<- lm(maturity.index ~ Longitude*orchard.type, data=c)
summary(mt3)
Anova(mt3)
#Longitude                6.807   1  4.5859 0.03433 *
#orchard.type             2.193   1  1.4777 0.22661  
#Longitude:orchard.type   0.157   1  0.1055 0.74589 


#Which abiotic factors are the most important drivers of fruit quality?----
####Looking at fixed variables here. Orchard num is nested within site code. Not focused on looking at mgmt systms 
#total precipitation x otype 
pl1<- glmmTMB(SSC ~ Szn.Total.Precip*orchard.type + (1|site.code/onum), data=c)
summary(pl1)
Anova(pl1)
#Szn.Total.Precip              7.3440  1   0.006729 **
#orchard.type                  1.5530  1   0.212698   
#Szn.Total.Precip:orchard.type 1.6851  1   0.194244  

pl2<- glmmTMB(Firmness ~ Szn.Total.Precip*orchard.type + (1|site.code/onum), data=c)
summary(pl2)
Anova(pl2)
#Szn.Total.Precip              0.1953  1     0.6586
#orchard.type                  0.0057  1     0.9397
#Szn.Total.Precip:orchard.type 0.6159  1     0.4326

pl3<- glmmTMB(avgwgt ~ Szn.Total.Precip*orchard.type + (1|site.code/onum), data=c)
summary(pl3)
Anova(pl3)
#Szn.Total.Precip              2.8708  1    0.09020 .
#orchard.type                  0.0067  1    0.93480  
#Szn.Total.Precip:orchard.type 3.9157  1    0.04784 *

pl4<- glmmTMB(maturity.index ~ Szn.Total.Precip*orchard.type + (1|site.code/onum), data=c)
summary(pl4)
Anova(pl4)
#Szn.Total.Precip              7.7214  1   0.005457 **
#orchard.type                  0.0916  1   0.762174   
#Szn.Total.Precip:orchard.type 0.2209  1   0.638342  

#avgerage seasonal temp x otype   
atl<- glmmTMB(SSC ~ Szn.Temp.Avg*orchard.type + (1|site.code/onum), data=c)
summary(atl)
Anova(atl)
#Szn.Temp.Avg              26.5282  1  2.597e-07 ***
#orchard.type               0.1087  1     0.7416    
#Szn.Temp.Avg:orchard.type  0.5756  1     0.4481  

at2<- glmmTMB(Firmness ~ Szn.Temp.Avg*orchard.type + (1|site.code/onum), data=c)
summary(at2)
Anova(at2)
#Szn.Temp.Avg              11.3947  1  0.0007365 ***
#orchard.type               0.0010  1  0.9748785    
#Szn.Temp.Avg:orchard.type  0.2247  1  0.6354858  

at3<- glmmTMB(avgwgt ~ Szn.Temp.Avg*orchard.type + (1|site.code/onum), data=c)
summary(at3)
Anova(at3)
#Szn.Temp.Avg              0.3241  1   0.569146   
#orchard.type              0.1482  1   0.700282   
#Szn.Temp.Avg:orchard.type 9.2367  1   0.002372 **

at4<- glmmTMB(maturity.index ~ Szn.Temp.Avg*orchard.type + (1|site.code/onum), data=c)
summary(at4)
Anova(at4)
#Szn.Temp.Avg              6.9050  1   0.008595 **
#orchard.type              0.0423  1   0.837056   
#Szn.Temp.Avg:orchard.type 0.6519  1   0.419438 

#proximity to water x otype 
pxl<- glmmTMB(SSC ~ Prox.Water*orchard.type + (1|site.code/onum), data=c)
summary(pxl)
Anova(pxl)
#Prox.Water              1.4996  1     0.2207
#orchard.type            0.9875  1     0.3204
#Prox.Water:orchard.type 0.4179  1     0.5180

px2<- glmmTMB(Firmness ~ Prox.Water*orchard.type + (1|site.code/onum), data=c)
summary(px2)
Anova(px2)
#Prox.Water              1.1828  1     0.2768
#orchard.type            0.0378  1     0.8457
#Prox.Water:orchard.type 1.8788  1     0.1705

px3<- glmmTMB(avgwgt ~ Prox.Water*orchard.type + (1|site.code/onum), data=c)
summary(px3)
Anova(px3)
#Prox.Water              0.1974  1     0.6569
#orchard.type            0.0106  1     0.9179
#Prox.Water:orchard.type 0.0027  1     0.9585

px4<- glmmTMB(maturity.index ~ Prox.Water*orchard.type + (1|site.code/onum), data=c)
summary(px4)
Anova(px4)
#NaN

#season avg max x otype 
smx1<- glmmTMB(SSC ~ Szn.Max.Avg*orchard.type + (1|site.code/onum), data=c)
summary(smx1)
Anova(smx1)
#Szn.Max.Avg              18.9030  1  1.375e-05 ***
#orchard.type              0.1568  1     0.6921    
#Szn.Max.Avg:orchard.type  0.0080  1     0.9287 

smx2<- glmmTMB(Firmness ~ orchard.type*Szn.Max.Avg + (1|site.code/onum), data=c)
summary(smx2)
Anova(smx2)
#orchard.type              0.0353  1   0.851057   
#Szn.Max.Avg              10.8113  1   0.001009 **
#orchard.type:Szn.Max.Avg  0.2284  1   0.632709   

smx3<- glmmTMB(avgwgt ~ orchard.type*Szn.Max.Avg + (1|site.code/onum), data=c)
summary(smx3)
Anova(smx3)
#orchard.type              0.1452  1  0.7031488    
#Szn.Max.Avg               0.0018  1  0.9663953    
#orchard.type:Szn.Max.Avg 11.4873  1  0.0007007 ***

smx4<- glmmTMB(maturity.index ~ orchard.type*Szn.Max.Avg + (1|site.code/onum), data=c)
summary(smx4)
Anova(smx4)
#orchard.type             0.3166  1   0.573643   
#Szn.Max.Avg              7.7421  1   0.005395 **
#orchard.type:Szn.Max.Avg 0.1453  1   0.703109 

#season avg min x otype 
sm1<- glmmTMB(SSC ~ orchard.type*Szn.Min.Avg + (1|site.code/onum), data=c)
summary(sm1)
Anova(sm1)
#orchard.type              1.0205  1     0.3124    
#Szn.Min.Avg              15.9763  1  6.414e-05 ***
#orchard.type:Szn.Min.Avg  0.4255  1     0.5142

sm2<- glmmTMB(Firmness ~ orchard.type*Szn.Min.Avg + (1|site.code/onum), data=c)
summary(sm2)
Anova(sm2)
#orchard.type             0.0786  1   0.779246   
#Szn.Min.Avg              9.3198  1   0.002267 **
#orchard.type:Szn.Min.Avg 0.0299  1   0.862609  

sm3<- glmmTMB(avgwgt ~ orchard.type*Szn.Min.Avg + (1|site.code/onum), data=c)
summary(sm3)
Anova(sm3)
#orchard.type             0.0461  1    0.82994  
#Szn.Min.Avg              0.8381  1    0.35995  
#orchard.type:Szn.Min.Avg 5.3243  1    0.02103 *

sm4<- glmmTMB(maturity.index ~ orchard.type*Szn.Min.Avg + (1|site.code/onum), data=c)
summary(sm4)
Anova(sm4)
#orchard.type             0.0811  1     0.7757
#Szn.Min.Avg              2.5495  1     0.1103
#orchard.type:Szn.Min.Avg 1.0085  1     0.3153


#Composite variable analysis----------------------------------------------------
library(car)
##looking at composite varibales 
c[12:16]
##plotting elationships 
scatterplotMatrix(c[12:16])





#Which specific management practices are the most important drivers of fruit chemistry and quality?------

names(c)

##pest pressure## 
pest1 <- glmmTMB(SSC ~ Pest_Index*orchard.type + (1|site.code/onum/Tree), data=c)
summary(pest1)
Anova(pest1)
#Pest_Index:orchard.type 2.9152  1    0.08775 .

pest2 <- glmmTMB(avgwgt ~ Pest_Index*orchard.type + (1|site.code/onum/Tree), data=c)
summary(pest2)
Anova(pest2)
# none

pest3 <- glmmTMB(Firmness ~ Pest_Index*orchard.type + (1|site.code/onum/Tree), data=c)
summary(pest3)
Anova(pest3)
#none 

pest4 <- glmmTMB(maturity.index ~ Pest_Index*orchard.type + (1|site.code/onum/Tree), data=c)
summary(pest4)
Anova(pest4)
#none 

##Herbicides## 
herb1 <- glmmTMB(SSC ~ Herbicides*orchard.type + (1|site.code/onum/Tree), data=c)
summary(herb1)
Anova(herb1)
#Herbicides              5.4980  1    0.01904 *
#orchard.type            2.7336  1    0.09826 .
#Herbicides:orchard.type 1.2220  1    0.26896  

herb2 <- glmmTMB(avgwgt ~ Herbicides*orchard.type + (1|site.code/onum/Tree), data=c)
summary(herb2)
Anova(herb2)
#none

herb3 <- glmmTMB(Firmness ~ Herbicides*orchard.type + (1|site.code/onum/Tree), data=c)
summary(herb3)
Anova(herb3)
#none 

herb4 <- glmmTMB(maturity.index ~ Herbicides*orchard.type + (1|site.code/onum/Tree), data=c)
summary(herb4)
Anova(herb4)
#none 


##Acres## 
Acres <- glmmTMB(SSC ~ Acres*orchard.type + (1|site.code/onum/Tree), data=c)
summary(Acres)
Anova(Acres)
#none

Acres <- glmmTMB(avgwgt ~ Acres*orchard.type + (1|site.code/onum/Tree), data=c)
summary(Acres)
Anova(Acres)
#Acres              10.8040  1   0.001013 **

Acres <- glmmTMB(Firmness ~ Acres*orchard.type + (1|site.code/onum/Tree), data=c)
summary(Acres)
Anova(Acres)
#none 

Acres <- glmmTMB(maturity.index ~ Acres*orchard.type + (1|site.code/onum/Tree), data=c)
summary(Acres)
Anova(Acres)
#none 

##Mowing## 
Mowing <- lm(SSC ~ Mowing*orchard.type, data=c)
summary(Mowing)
Anova(Mowing)
#none

Mowing <- lm(avgwgt ~ Mowing*orchard.type, data=c)
summary(Mowing)
Anova(Mowing)
#none 

Mowing <- lm(Firmness ~ Mowing*orchard.type, data=c)
summary(Mowing)
Anova(Mowing)
#none 

Mowing <- lm(maturity.index ~ Mowing*orchard.type, data=c)
summary(Mowing)
Anova(Mowing)
#none 





